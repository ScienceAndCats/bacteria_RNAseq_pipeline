#!/usr/bin/env python3
"""Submit isolated per-sample Slurm jobs and publish combined project reports."""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import uuid

from sample_manifest import build_manifest
from merge_featurecounts import merge_featurecounts

REPOSITORY = Path(__file__).resolve().parent
SETTINGS = (
    'FASTQ_DIR', 'FASTQ_GLOB', 'READ_LAYOUT', 'OUTPUT_DIR', 'MAX_SAMPLES',
    'ARRAY_CONCURRENCY', 'THREADS', 'REFERENCE_CACHE_DIR', 'CONDA_SH', 'CONDA_ENV',
    'SLURM_ACCOUNT', 'SLURM_MAIL_USER',
    'DECOY_REFERENCE_DIR', 'BACTERIA_REFERENCE_DIR', 'HOST_REFERENCE_DIR',
    'ADAPTER_SINGLE', 'ADAPTER_R1', 'ADAPTER_R2', 'ADAPTER_SEQUENCE',
    'MIN_READ_LENGTH', 'SAMPLE_READS', 'SAMPLE_SIZE', 'SAMPLE_SEED',
    'GENE_POSITION_BINS', 'METAGENE_MIN_FEATURE_READS', 'COUNT_STRAND',
    'CSV_CONVERSION_SCRIPT',
)
DEFAULTS = {
    'FASTQ_DIR': '.', 'FASTQ_GLOB': '*.fastq*', 'READ_LAYOUT': 'single',
    'OUTPUT_DIR': '.', 'MAX_SAMPLES': '40', 'ARRAY_CONCURRENCY': '8',
    'SLURM_ACCOUNT': '', 'SLURM_MAIL_USER': '',
    'THREADS': '16', 'REFERENCE_CACHE_DIR': '', 'CONDA_SH': '',
    'CONDA_ENV': 'bacteria-RNAseq-pipeline', 'HOST_REFERENCE_DIR': '',
    'CSV_CONVERSION_SCRIPT': 'bacteria_with_decoys_csvConversion.py',
}


def positive_integer(value, name):
    if not re.fullmatch(r'[1-9][0-9]*', str(value)):
        raise ValueError(f'{name} must be a positive integer')
    return int(value)


def absolute_path(value, base):
    if any(ord(char) < 32 or ord(char) == 127 for char in str(value)):
        raise ValueError('Paths must not contain control characters')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def write_json(path, value):
    path = Path(path)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, prefix=f'.{path.name}.', delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2)
        handle.write('\n')
    temporary.replace(path)


def write_config(path, settings):
    path.write_text('# Frozen settings generated at submission; edit the source config for a new run.\n' + ''.join(f'{key}={shlex.quote(str(value))}\n' for key, value in settings.items()))


def create_plan(config, base_dir, environment=None):
    environment = os.environ if environment is None else environment
    base = Path(base_dir).resolve()
    settings = {key: environment[key] for key in SETTINGS if key in environment}
    for key, value in DEFAULTS.items():
        settings.setdefault(key, value)
    maximum = positive_integer(settings['MAX_SAMPLES'], 'MAX_SAMPLES')
    concurrency = positive_integer(settings['ARRAY_CONCURRENCY'], 'ARRAY_CONCURRENCY')
    threads = positive_integer(settings['THREADS'], 'THREADS')
    for key in ('DECOY_REFERENCE_DIR', 'BACTERIA_REFERENCE_DIR'):
        if not settings.get(key):
            raise ValueError(f'Set {key} to the reference folder in the config')
    settings['FASTQ_DIR'] = str(absolute_path(settings['FASTQ_DIR'], base))
    settings['OUTPUT_DIR'] = str(absolute_path(settings['OUTPUT_DIR'], base))
    for key in ('DECOY_REFERENCE_DIR', 'BACTERIA_REFERENCE_DIR', 'HOST_REFERENCE_DIR', 'CONDA_SH'):
        if settings.get(key):
            settings[key] = str(absolute_path(settings[key], base))
    if settings['CONDA_SH'] and not Path(settings['CONDA_SH']).is_file():
        raise ValueError(f'CONDA_SH does not exist: {settings["CONDA_SH"]}')
    settings['REFERENCE_CACHE_DIR'] = str(absolute_path(settings['REFERENCE_CACHE_DIR'], base)) if settings['REFERENCE_CACHE_DIR'] else str(Path(settings['OUTPUT_DIR']) / '.reference_cache')
    settings['CSV_CONVERSION_SCRIPT'] = str(absolute_path(settings['CSV_CONVERSION_SCRIPT'], REPOSITORY))
    if not Path(settings['CSV_CONVERSION_SCRIPT']).is_file():
        raise ValueError(f'CSV_CONVERSION_SCRIPT does not exist: {settings["CSV_CONVERSION_SCRIPT"]}')
    manifest = build_manifest(settings['FASTQ_DIR'], settings['FASTQ_GLOB'], settings['READ_LAYOUT'], max_samples=maximum)
    output = Path(settings['OUTPUT_DIR'])
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:10]
    run_dir = output / 'slurm_runs' / run_id
    run_dir.mkdir(parents=True)
    (run_dir / 'logs').mkdir()
    sample_manifest = run_dir / 'samples.json'
    write_json(sample_manifest, manifest)
    frozen = run_dir / 'run.env'
    write_config(frozen, settings)
    prepare_settings = dict(settings, OUTPUT_DIR=str(run_dir))
    prepare_config = run_dir / 'prepare.env'
    write_config(prepare_config, prepare_settings)
    samples = []
    for index, sample in enumerate(manifest['samples']):
        sample_dir = run_dir / 'samples' / f'{index:04d}_{sample["sample"]}'
        sample_dir.mkdir(parents=True)
        sample_config = sample_dir / 'run.env'
        write_config(sample_config, dict(settings, OUTPUT_DIR=str(sample_dir)))
        samples.append(dict(sample, output_dir=str(sample_dir), config=str(sample_config)))
    template = run_dir / 'run_bacteria_pipeline.sbatch'
    shutil.copyfile(REPOSITORY / 'run_bacteria_pipeline.sbatch', template)
    plan = {
        'version': 1, 'run_id': run_id, 'run_dir': str(run_dir),
        'output_dir': str(output), 'base_dir': str(base),
        'source_config': str(absolute_path(config, base)), 'frozen_config': str(frozen),
        'prepare_config': str(prepare_config), 'sample_manifest': str(sample_manifest),
        'reference_manifest': str(run_dir / 'reference_manifest.json'),
        'read_layout': manifest['read_layout'], 'samples': samples,
        'threads': threads, 'array_concurrency': concurrency, 'max_samples': maximum,
        'slurm_account': settings['SLURM_ACCOUNT'],
        'slurm_mail_user': settings['SLURM_MAIL_USER'],
        'pipeline_script': str(REPOSITORY / 'map_bacteria_with_decoys.sh'),
        'workflow_script': str(REPOSITORY / 'slurm_workflow.py'),
        'sbatch_script': str(template), 'csv_script': settings['CSV_CONVERSION_SCRIPT'],
    }
    path = run_dir / 'plan.json'
    write_json(path, plan)
    return path, plan


def load_plan(path):
    plan = json.loads(Path(path).read_text())
    if plan.get('version') != 1 or not plan.get('samples'):
        raise ValueError('Unsupported or empty Slurm run plan')
    return plan


def scheduler_command(plan, path, mode, dependency=None):
    log_directory = str(Path(plan['run_dir']) / 'logs').replace('%', '%%')
    command = ['sbatch', '--parsable', f'--job-name=bacteria_{mode}', f'--chdir={plan["base_dir"]}', f'--output={log_directory}/%x_%A_%a.out']
    if plan.get('slurm_account'):
        command.append(f'--account={plan["slurm_account"]}')
    if plan.get('slurm_mail_user'):
        command += [f'--mail-user={plan["slurm_mail_user"]}', '--mail-type=BEGIN,END,FAIL']
    else:
        command.append('--mail-type=NONE')
    if mode == 'finalize':
        command += ['--cpus-per-task=1', '--mem=4G']
    else:
        command += [f'--cpus-per-task={plan["threads"]}']
    if mode == 'sample':
        command += [f'--array=0-{len(plan["samples"]) - 1}%{plan["array_concurrency"]}']
    if dependency:
        command += [f'--dependency=afterok:{dependency}', '--kill-on-invalid-dep=yes']
    return command + [plan['sbatch_script'], mode, str(path), plan['workflow_script'], plan['frozen_config']]


def submit_plan(path, plan):
    ids = {}
    try:
        for mode, predecessor in (('prepare', None), ('sample', 'prepare'), ('finalize', 'sample')):
            command = scheduler_command(plan, path, mode, ids.get(predecessor))
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            response = result.stdout.strip()
            if not re.fullmatch(r'[0-9]+(?:;[A-Za-z0-9_.-]+)?', response):
                raise ValueError(f'sbatch returned an unexpected job ID for {mode}: {response!r}')
            ids[mode] = response.split(';', 1)[0]
            write_json(Path(plan['run_dir']) / 'submission.json', {'status': 'submitting', 'jobs': ids})
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        cancelled = False
        if ids and shutil.which('scancel'):
            cancelled = subprocess.run(['scancel', *ids.values()], capture_output=True, text=True).returncode == 0
        write_json(Path(plan['run_dir']) / 'submission.json', {'status': 'failed', 'jobs': ids, 'cancellation_requested': cancelled})
        detail = exc.stderr.strip() if isinstance(exc, subprocess.CalledProcessError) and exc.stderr else str(exc)
        state = 'Cancellation requested for' if cancelled else 'Check/cancel these already submitted jobs:'
        raise ValueError(f'Job submission failed: {detail}. {state} {ids or "none"}') from exc
    write_json(Path(plan['run_dir']) / 'submission.json', {'status': 'submitted', 'jobs': ids})
    return ids


def expected_stages(references):
    names = ['decoy']
    for kind in ('rna', 'dna'):
        for organism in ('bacteria', 'host'):
            if references.get(organism) and references[organism].get(kind):
                names.append(f'{organism}_{kind}')
    return names


def required_file(value, directory, label):
    if not value:
        raise ValueError(f'Missing {label} path')
    path = Path(value)
    path = path if path.is_absolute() else directory / path
    path = path.resolve()
    if not path.is_file():
        raise ValueError(f'Missing {label}: {path}')
    return str(path)


def collect_samples(plan):
    references = json.loads(Path(plan['reference_manifest']).read_text())
    expected = expected_stages(references)
    combined = []
    for index, sample in enumerate(plan['samples']):
        directory = Path(sample['output_dir'])
        worker = json.loads((directory / 'run_manifest.json').read_text())
        if worker.get('read_layout') != plan['read_layout']:
            raise ValueError(f'Read layout mismatch for {sample["sample"]}')
        if worker.get('sample_index') != index or worker.get('input_sample_manifest') != plan['sample_manifest']:
            raise ValueError(f'Array manifest provenance mismatch for {sample["sample"]}')
        entries = worker.get('samples', [])
        if len(entries) != 1 or entries[0].get('sample') != sample['sample']:
            raise ValueError(f'Missing, duplicate, or unexpected sample in worker {index}')
        entry = dict(entries[0])
        entry['cutadapt_log'] = required_file(entry.get('cutadapt_log'), directory, 'cutadapt log')
        stages = [dict(stage) for stage in entry.get('stages', [])]
        if [stage.get('name') for stage in stages] != expected:
            raise ValueError(f'Incomplete or unexpected stages for {sample["sample"]}')
        for stage in stages:
            stage['log'] = required_file(stage.get('log'), directory, 'Bowtie2 log')
            if stage.get('organism') != 'decoy':
                stage['categories'] = required_file(stage.get('categories'), directory, 'category report')
                stage['coverage'] = required_file(stage.get('coverage'), directory, 'coverage report')
        entry['stages'] = stages
        leftovers = entry.get('leftover', [])
        if len(leftovers) != (2 if plan['read_layout'] == 'paired' else 1):
            raise ValueError(f'Incomplete leftover FASTQs for {sample["sample"]}')
        entry['leftover'] = [required_file(path, directory, 'leftover FASTQ') for path in leftovers]
        combined.append(entry)
    return {'read_layout': plan['read_layout'], 'slurm_run_id': plan['run_id'], 'samples': combined}


def validate_gene_tables(plan):
    """Require all featureCounts tables implied by prepared annotations.

    The merge utility also accepts workflows with no given feature type. Here
    the reference manifest supplies the expected stage/type matrix, so losing a
    type from every worker cannot silently produce an incomplete project.
    """
    references_path = Path(plan['reference_manifest'])
    references = json.loads(references_path.read_text())
    expected = []
    for organism in ('bacteria', 'host'):
        for kind in ('rna', 'dna'):
            stage = (references.get(organism) or {}).get(kind)
            if not stage:
                continue
            gff = Path(required_file(stage.get('gff'), references_path.parent, 'prepared annotation'))
            with gff.open() as annotations:
                types = {line.split('\t')[2] for line in annotations if line.strip() and not line.startswith('#')}
            directory = Path('.') if organism == 'bacteria' and kind == 'dna' else Path('bowtie_alignments') / organism
            if kind == 'rna':
                directory /= 'rna'
            for feature in ('CDS', 'rRNA'):
                if feature not in types:
                    continue
                suffix = 'summary' if feature == 'CDS' else 'rRNA'
                table = directory / f'featurecounts_{organism.upper()}_{suffix}.txt'
                expected.extend((table, Path(str(table) + '.summary')))
    for sample in plan['samples']:
        directory = Path(sample['output_dir'])
        for relative in expected:
            required_file(str(relative), directory, f'featureCounts output for {sample["sample"]}')



@contextmanager
def publication_lock(output):
    with (output / '.slurm_finalize.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def finalize(plan):
    output = Path(plan['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    with publication_lock(output):
        manifest = collect_samples(plan)
        validate_gene_tables(plan)
        staging = Path(tempfile.mkdtemp(prefix='final_reports-', dir=plan['run_dir']))
        try:
            write_json(staging / 'run_manifest.json', manifest)
            summary_name = f'BACTERIA_{output.name}.csv'
            alias_name = f'BACTERIA_{output.name}_read_breakdown.csv'
            subprocess.run([sys.executable, plan['csv_script'], '--output-dir', str(staging), '--output', str(staging / summary_name)], check=True)
            subprocess.run([sys.executable, str(REPOSITORY / 'read_breakdown_csv.py'), '--output-dir', str(staging), '--output', str(staging / 'output_read_breakdown.csv')], check=True)
            shutil.copyfile(staging / 'output_read_breakdown.csv', staging / alias_name)
            shutil.copyfile(plan['reference_manifest'], staging / 'reference_manifest.json')
            gene_tables = merge_featurecounts([Path(sample['output_dir']) for sample in plan['samples']], staging)
            names = [summary_name, 'output_read_breakdown.csv', alias_name, 'reference_manifest.json']
            names += [str(path.relative_to(staging)) for path in gene_tables]
            names.append('run_manifest.json')
            published = output / '.slurm_published_files.json'
            prior = json.loads(published.read_text()).get('files', []) if published.exists() else []
            if not isinstance(prior, list) or any(not isinstance(name, str) or Path(name).is_absolute() or '..' in Path(name).parts for name in prior):
                raise ValueError('Invalid previous Slurm publication metadata')
            # All validation/report generation finishes before replacing the last
            # successful project's reports. Unique staging prevents stale rows.
            for name in names:
                destination = output / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = destination.with_name(f'.{destination.name}.{plan["run_id"]}.tmp')
                shutil.copyfile(staging / name, temporary)
                temporary.replace(destination)
            # Remove only files tracked as generated by a previous successful
            # finalizer, e.g. host tables after a new run disables the host.
            for name in set(prior) - set(names):
                (output / name).unlink(missing_ok=True)
            write_json(published, {'version': 1, 'run_id': plan['run_id'], 'files': names})
            write_json(Path(plan['run_dir']) / 'finalized.json', {'status': 'completed', 'samples': len(manifest['samples']), 'output_dir': str(output)})
        finally:
            shutil.rmtree(staging)


def execute(mode, plan):
    if mode == 'finalize':
        finalize(plan)
        return
    command = ['bash', plan['pipeline_script']]
    if mode == 'prepare':
        command += [plan['prepare_config'], '--prepare-only']
    elif mode == 'sample':
        value = os.environ.get('SLURM_ARRAY_TASK_ID', '')
        if not re.fullmatch(r'[0-9]+', value) or not 0 <= int(value) < len(plan['samples']):
            raise ValueError('SLURM_ARRAY_TASK_ID must identify a sample in the frozen manifest')
        sample = plan['samples'][int(value)]
        (Path(sample['output_dir']) / 'run_manifest.json').unlink(missing_ok=True)
        command += [sample['config'], '--sample-manifest', plan['sample_manifest'], '--sample-index', value, '--prepared-references', plan['reference_manifest'], '--skip-reports']
    else:
        raise ValueError(f'Unknown job mode: {mode}')
    subprocess.run(command, cwd=plan['base_dir'], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    submit = commands.add_parser('submit', help='Freeze samples/settings and submit the dependency chain')
    submit.add_argument('--config', type=Path, required=True)
    submit.add_argument('--base-dir', type=Path, default=Path.cwd())
    submit.add_argument('--dry-run', action='store_true')
    configuration = commands.add_parser('configuration', help='Print frozen config path for sbatch activation')
    configuration.add_argument('plan', type=Path)
    run = commands.add_parser('execute', help='Run one scheduled stage')
    run.add_argument('mode', choices=('prepare', 'sample', 'finalize'))
    run.add_argument('plan', type=Path)
    final = commands.add_parser('finalize', help='Validate all workers and regenerate combined CSVs')
    final.add_argument('plan', type=Path)
    args = parser.parse_args()
    try:
        if args.command == 'submit':
            path, plan = create_plan(args.config, args.base_dir)
            print(f'Frozen plan: {path}', flush=True)
            print(f'{len(plan["samples"])} samples; array 0-{len(plan["samples"]) - 1}%{plan["array_concurrency"]}', flush=True)
            if args.dry_run:
                for mode, dependency in (('prepare', None), ('sample', '<prepare-job-id>'), ('finalize', '<array-job-id>')):
                    print(shlex.join(scheduler_command(plan, path, mode, dependency)))
            else:
                ids = submit_plan(path, plan)
                print(f'Submitted prepare={ids["prepare"]}, array={ids["sample"]}, final={ids["finalize"]}')
        else:
            plan = load_plan(args.plan)
            if args.command == 'configuration':
                print(plan['frozen_config'])
            elif args.command == 'execute':
                execute(args.mode, plan)
            else:
                finalize(plan)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
