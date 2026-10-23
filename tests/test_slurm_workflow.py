"""Slurm planning, submission dependencies, and finalization safety checks."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import slurm_workflow as workflow


class SlurmWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='slurm workflow ')
        self.base = Path(self.temp.name)
        self.reads = self.base / 'read files'
        self.reads.mkdir()
        self.output = self.base / 'analysis results'
        self.config = self.base / 'original.env'
        self.config.write_text('# source settings\n')
        for name in ('zebra', 'alpha'):
            (self.reads / f'{name}.fastq').write_text('@read\nACGT\n+\nIIII\n')
        self.environment = {
            'FASTQ_DIR': str(self.reads), 'OUTPUT_DIR': str(self.output),
            'DECOY_REFERENCE_DIR': str(self.base / 'decoy'),
            'BACTERIA_REFERENCE_DIR': str(self.base / 'bacteria'),
        }

    def tearDown(self):
        self.temp.cleanup()

    def plan(self, **overrides):
        return workflow.create_plan(self.config, self.base, dict(self.environment, **overrides))

    def test_freezes_sorted_samples_resolved_settings_and_unique_runs(self):
        path, plan = self.plan(THREADS='3', ADAPTER_SEQUENCE='legacy adapter')
        frozen = Path(plan['frozen_config']).read_text()
        self.assertIn('ADAPTER_SEQUENCE=', frozen)
        self.assertIn('legacy adapter', frozen)
        self.assertIn('THREADS=3', frozen)
        self.assertEqual([s['sample'] for s in plan['samples']], ['alpha', 'zebra'])
        self.assertTrue(Path(plan['samples'][0]['output_dir']).is_dir())
        self.assertNotEqual(plan['samples'][0]['output_dir'], plan['samples'][1]['output_dir'])
        self.assertEqual(plan['max_samples'], 40)
        self.assertEqual(plan['array_concurrency'], 8)
        self.assertEqual(json.loads(path.read_text()), plan)
        self.config.write_text('OUTPUT_DIR=changed\n')
        (self.reads / 'new.fastq').write_text('')
        self.assertEqual(len(json.loads(Path(plan['sample_manifest']).read_text())['samples']), 2)
        _, second = self.plan()
        self.assertNotEqual(plan['run_dir'], second['run_dir'])
        self.assertEqual(frozen, Path(plan['frozen_config']).read_text())

    def test_default_rejects_more_than_40_without_silent_truncation(self):
        for index in range(39):
            (self.reads / f'extra_{index}.fastq').write_text('')
        with self.assertRaises(ValueError):
            self.plan()
        self.assertFalse(self.output.exists())
        _, plan = self.plan(MAX_SAMPLES='41')
        self.assertEqual(len(plan['samples']), 41)

    def test_rejects_invalid_limits_and_path_controls(self):
        for changes in ({'MAX_SAMPLES': '0'}, {'ARRAY_CONCURRENCY': '-1'}, {'THREADS': 'nan'}, {'OUTPUT_DIR': 'bad\nname'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.plan(**changes)

    def test_submit_dependencies_and_exact_array_size(self):
        path, plan = self.plan(SLURM_ACCOUNT='research_lab', SLURM_MAIL_USER='researcher@example.org')
        calls = []
        def fake(command, **kwargs):
            calls.append(command)
            return subprocess.CompletedProcess(command, 0, str(100 + len(calls)) + ';cluster\n', '')
        with patch.object(workflow.subprocess, 'run', side_effect=fake):
            ids = workflow.submit_plan(path, plan)
        self.assertEqual(ids, {'prepare': '101', 'sample': '102', 'finalize': '103'})
        self.assertFalse(any('--array' in value for value in calls[0]))
        self.assertIn('--array=0-1%8', calls[1])
        self.assertIn('--dependency=afterok:101', calls[1])
        self.assertIn('--dependency=afterok:102', calls[2])
        self.assertIn('--kill-on-invalid-dep=yes', calls[1])
        self.assertIn('--kill-on-invalid-dep=yes', calls[2])
        for command in calls:
            self.assertIn('--account=research_lab', command)
            self.assertIn('--mail-user=researcher@example.org', command)
            self.assertIn('--mail-type=BEGIN,END,FAIL', command)
        for config in (plan['frozen_config'], plan['prepare_config'], *(sample['config'] for sample in plan['samples'])):
            frozen = Path(config).read_text()
            self.assertIn('SLURM_ACCOUNT=research_lab\n', frozen)
            self.assertIn('SLURM_MAIL_USER=researcher@example.org\n', frozen)
        self.assertEqual(calls[0][-1], plan['frozen_config'])
        self.assertEqual(json.loads((Path(plan['run_dir']) / 'submission.json').read_text())['status'], 'submitted')

    def test_empty_account_and_email_use_defaults_without_notifications(self):
        path, plan = self.plan()
        for mode in ('prepare', 'sample', 'finalize'):
            command = workflow.scheduler_command(plan, path, mode)
            self.assertFalse(any(value.startswith('--account=') for value in command))
            self.assertFalse(any(value.startswith('--mail-user=') for value in command))
            self.assertIn('--mail-type=NONE', command)

    def test_account_and_email_can_be_configured_independently(self):
        for settings in ({'SLURM_ACCOUNT': 'research_lab'}, {'SLURM_MAIL_USER': 'researcher@example.org'}):
            with self.subTest(settings=settings):
                path, plan = self.plan(**settings)
                command = workflow.scheduler_command(plan, path, 'sample')
                self.assertEqual('--account=research_lab' in command, 'SLURM_ACCOUNT' in settings)
                self.assertEqual('--mail-user=researcher@example.org' in command, 'SLURM_MAIL_USER' in settings)
                self.assertEqual('--mail-type=BEGIN,END,FAIL' in command, 'SLURM_MAIL_USER' in settings)

    def test_partial_submission_failure_cancels_only_own_jobs(self):
        path, plan = self.plan()
        calls = []
        def fake(command, **kwargs):
            calls.append(command)
            if command[0] == 'scancel':
                return subprocess.CompletedProcess(command, 0, '', '')
            if len(calls) == 1:
                return subprocess.CompletedProcess(command, 0, '715\n', '')
            raise subprocess.CalledProcessError(1, command, stderr='account rejected')
        with patch.object(workflow.subprocess, 'run', side_effect=fake), patch.object(workflow.shutil, 'which', return_value='/bin/scancel'):
            with self.assertRaisesRegex(ValueError, 'account rejected'):
                workflow.submit_plan(path, plan)
        self.assertEqual(calls[-1], ['scancel', '715'])
        state = json.loads((Path(plan['run_dir']) / 'submission.json').read_text())
        self.assertEqual(state['jobs'], {'prepare': '715'})
        self.assertTrue(state['cancellation_requested'])

    def test_sample_execute_uses_exact_index_and_invalidates_old_manifest(self):
        _, plan = self.plan()
        old = Path(plan['samples'][1]['output_dir']) / 'run_manifest.json'
        old.write_text('stale')
        with patch.dict(os.environ, {'SLURM_ARRAY_TASK_ID': '1'}), patch.object(workflow.subprocess, 'run') as run:
            workflow.execute('sample', plan)
        self.assertFalse(old.exists())
        command = run.call_args.args[0]
        self.assertEqual(command[2], plan['samples'][1]['config'])
        self.assertIn('--skip-reports', command)
        self.assertEqual(command[command.index('--sample-index') + 1], '1')
        for value in ('-1', '2', 'nope', ''):
            with patch.dict(os.environ, {'SLURM_ARRAY_TASK_ID': value}), self.assertRaises(ValueError):
                workflow.execute('sample', plan)

    def test_log_path_percent_is_escaped_but_array_placeholders_are_not(self):
        path, plan = self.plan()
        plan['run_dir'] += '%j'
        command = workflow.scheduler_command(plan, path, 'sample')
        log = next(value for value in command if value.startswith('--output='))
        self.assertIn('%%j/logs/%x_%A_%a.out', log)

    def test_submit_wrapper_dry_run_uses_config_folder_and_submits_nothing(self):
        self.config.write_text("FASTQ_DIR='read files'\nOUTPUT_DIR='analysis results'\nDECOY_REFERENCE_DIR=decoy\nBACTERIA_REFERENCE_DIR=bacteria\nTHREADS=2\nSLURM_ACCOUNT=research_lab\nSLURM_MAIL_USER=researcher@example.org\n")
        result = subprocess.run(['bash', str(workflow.REPOSITORY / 'submit_bacteria_pipeline.sh'), str(self.config), '--dry-run'], cwd=self.base, text=True, capture_output=True, check=True)
        self.assertIn('2 samples; array 0-1%8', result.stdout)
        plans = list(self.output.glob('slurm_runs/*/plan.json'))
        self.assertEqual(len(plans), 1)
        plan = json.loads(plans[0].read_text())
        self.assertEqual(plan['samples'][0]['read1'], str(self.reads / 'alpha.fastq'))
        self.assertEqual(plan['slurm_account'], 'research_lab')
        self.assertEqual(plan['slurm_mail_user'], 'researcher@example.org')
        commands = [line for line in result.stdout.splitlines() if line.startswith('sbatch ')]
        self.assertEqual(len(commands), 3)
        for command in commands:
            self.assertIn('--account=research_lab', command)
            self.assertIn('--mail-user=researcher@example.org', command)
            self.assertIn('--mail-type=BEGIN,END,FAIL', command)
        self.assertFalse((Path(plan['run_dir']) / 'submission.json').exists())

    def test_sbatch_activation_runs_before_python_is_needed(self):
        activation = self.base / 'conda setup.sh'
        activation.write_text('conda() { export PATH=' + str(Path(os.sys.executable).parent) + ':/usr/bin:/bin; }\n')
        frozen = self.base / 'frozen.env'
        workflow.write_config(frozen, {'CONDA_SH': str(activation), 'CONDA_ENV': 'test-env'})
        stub = self.base / 'stub.py'
        stub.write_text("import sys; print('dispatch:' + sys.argv[1] + ':' + sys.argv[2])\n")
        result = subprocess.run(['/bin/bash', str(workflow.REPOSITORY / 'run_bacteria_pipeline.sbatch'), 'prepare', str(self.base / 'plan.json'), str(stub), str(frozen)], env={'PATH': ''}, capture_output=True, text=True, check=True)
        self.assertIn('dispatch:execute:prepare', result.stdout)

    def worker_manifests(self, plan):
        Path(plan['reference_manifest']).write_text(json.dumps({'decoy': {'index': 'decoy'}, 'bacteria': {'dna': {'index': 'dna', 'gff': 'genes'}}, 'host': None}))
        for index, sample in enumerate(plan['samples']):
            directory = Path(sample['output_dir'])
            files = {}
            for name in ('trim.log', 'decoy.log', 'dna.log', 'categories.json', 'coverage.tsv', 'leftover.fastq.gz'):
                path = directory / name
                path.write_text('report')
                files[name] = str(path)
            worker = {
                'read_layout': plan['read_layout'], 'sample_index': index,
                'input_sample_manifest': plan['sample_manifest'], 'samples': [{
                    'sample': sample['sample'], 'cutadapt_log': files['trim.log'],
                    'stages': [{'name': 'decoy', 'organism': 'decoy', 'kind': 'decoy', 'log': files['decoy.log']}, {'name': 'bacteria_dna', 'organism': 'bacteria', 'kind': 'dna', 'log': files['dna.log'], 'categories': files['categories.json'], 'coverage': files['coverage.tsv']}],
                    'leftover': [files['leftover.fastq.gz']],
                }],
            }
            (directory / 'run_manifest.json').write_text(json.dumps(worker))

    def test_collect_requires_all_workers_and_matching_frozen_provenance(self):
        _, plan = self.plan()
        self.worker_manifests(plan)
        combined = workflow.collect_samples(plan)
        self.assertEqual([sample['sample'] for sample in combined['samples']], ['alpha', 'zebra'])
        manifest = Path(plan['samples'][1]['output_dir']) / 'run_manifest.json'
        contents = json.loads(manifest.read_text())
        contents['sample_index'] = 0
        manifest.write_text(json.dumps(contents))
        with self.assertRaisesRegex(ValueError, 'provenance'):
            workflow.collect_samples(plan)
        manifest.unlink()
        with self.assertRaises(FileNotFoundError):
            workflow.collect_samples(plan)

    def test_missing_required_cds_in_every_worker_blocks_final_reports(self):
        _, plan = self.plan()
        self.worker_manifests(plan)
        gff = Path(plan['run_dir']) / 'genes.gff'
        gff.write_text('##gff-version 3\nchr\ttest\tCDS\t1\t100\t.\t+\t.\tlocus_tag=g1\n')
        references = json.loads(Path(plan['reference_manifest']).read_text())
        references['bacteria']['dna']['gff'] = str(gff)
        Path(plan['reference_manifest']).write_text(json.dumps(references))
        old = self.output / 'output_read_breakdown.csv'
        old.write_text('prior successful run\n')
        with self.assertRaisesRegex(ValueError, 'featureCounts output'):
            workflow.finalize(plan)
        self.assertEqual(old.read_text(), 'prior successful run\n')
        for sample in plan['samples']:
            table = Path(sample['output_dir']) / 'featurecounts_BACTERIA_summary.txt'
            table.write_text('generated table')
            Path(str(table) + '.summary').write_text('generated summary')
        workflow.validate_gene_tables(plan)
        Path(str(table) + '.summary').unlink()
        with self.assertRaisesRegex(ValueError, 'featureCounts output'):
            workflow.validate_gene_tables(plan)

    def test_failed_finalization_preserves_last_successful_reports(self):
        _, plan = self.plan()
        old = self.output / 'output_read_breakdown.csv'
        old.write_text('prior successful run\n')
        with self.assertRaises(FileNotFoundError):
            workflow.finalize(plan)
        self.assertEqual(old.read_text(), 'prior successful run\n')


if __name__ == '__main__':
    unittest.main()
