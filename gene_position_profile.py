#!/usr/bin/env python3
"""Create strand-aware feature-relative bacterial read/fragment profiles."""
import argparse, csv, gzip, os, statistics, sys
from collections import Counter, defaultdict
from dataclasses import dataclass

try:
    import pysam
except ImportError as exc:
    raise SystemExit("ERROR: pysam is required (install the pipeline conda environment).") from exc

@dataclass(frozen=True)
class Feature:
    id: str; name: str; type: str; cls: str; ref: str; start: int; end: int; strand: str

def attrs(text):
    out = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item: continue
        if "=" in item: key, value = item.split("=", 1)
        elif " " in item: key, value = item.split(None, 1); value = value.strip('"')
        else: continue
        out[key] = value.strip().strip('"')
    return out

def first(a, keys):
    return next((a[k].split(",")[0] for k in keys if a.get(k)), "")

def parse_gff(path, preferred):
    opener = gzip.open if path.endswith(".gz") else open
    records, by_id = [], {}
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"): continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9: continue
            try: start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError: raise SystemExit(f"ERROR: invalid GFF coordinates at {path}:{number}")
            if fields[6] not in ("+", "-") or end <= start: continue
            a = attrs(fields[8]); rec = {"ref":fields[0], "type":fields[2], "start":start,
                "end":end, "strand":fields[6], "a":a, "line":number}
            records.append(rec)
            if a.get("ID"): by_id[a["ID"]] = rec
    def ancestors(rec):
        seen=set(); stack=rec["a"].get("Parent", "").split(",") if rec["a"].get("Parent") else []
        while stack:
            ident=stack.pop()
            if ident in seen: continue
            seen.add(ident); parent=by_id.get(ident)
            if parent:
                yield parent
                stack.extend(parent["a"].get("Parent", "").split(",") if parent["a"].get("Parent") else [])
    features=[]; used=set()
    genes_by_key = {}
    for rec in records:
        if rec["type"].lower() == "gene":
            for key_name in (preferred, "locus", "locus_tag", "gene", "ID", "Name"):
                if rec["a"].get(key_name):
                    value = rec["a"][key_name].split(",")[0]
                    genes_by_key[(rec["ref"], value)] = rec
    cds_groups=defaultdict(list)
    for rec in records:
        if rec["type"].lower() != "cds": continue
        chain=list(ancestors(rec)); gene=next((x for x in chain if x["type"].lower()=="gene"), None)
        matching_gene = None
        if gene is None:
            for value in (first(rec["a"], [preferred]), first(rec["a"], ["locus", "locus_tag", "gene"])):
                if value and (rec["ref"], value) in genes_by_key:
                    matching_gene = genes_by_key[(rec["ref"], value)]; break
        owner=gene or matching_gene or (chain[-1] if chain else rec)
        key=first(owner["a"], [preferred,"locus","locus_tag","gene","ID","Name"]) or first(rec["a"], [preferred,"locus","locus_tag","gene","Parent","ID","Name"])
        if key: cds_groups[(owner["ref"], key)].append((owner, rec))
    for (_, key), group in cds_groups.items():
        owner=group[0][0]; cds=[x[1] for x in group]
        # A gene record supplies the complete biological gene interval; otherwise merge its CDS pieces.
        is_gene=owner["type"].lower()=="gene"
        start=owner["start"] if is_gene else min(x["start"] for x in cds)
        end=owner["end"] if is_gene else max(x["end"] for x in cds)
        a=owner["a"]; name=first(a,["gene","Name",preferred,"locus_tag","locus","ID"]) or key
        f=Feature(key,name,"gene" if is_gene else "CDS","CDS",owner["ref"],start,end,owner["strand"])
        features.append(f); used.add((f.cls,f.ref,f.start,f.end,f.strand))
    known={"rrna","trna","tmrna","ncrna","srna","misc_rna","snrna","scrna","snorna","mirna","pirna","antisense_rna","rnase_p_rna","rnase_mrp_rna","telomerase_rna","guide_rna","ribozyme"}
    for rec in records:
        low=rec["type"].lower()
        # RNA children are biological non-CDS units. Exclude generic genes and protein mRNA/transcripts.
        is_nc = low in known or ("rna" in low and low not in {"mrna","transcript","primary_transcript"})
        if not is_nc: continue
        a=rec["a"]; ident=first(a,[preferred,"locus","locus_tag","gene","ID","Name","Parent"])
        if not ident: continue
        key=("non_CDS",rec["ref"],rec["start"],rec["end"],rec["strand"])
        if key in used: continue
        name=first(a,["gene","Name",preferred,"locus_tag","locus","ID"]) or ident
        features.append(Feature(ident,name,rec["type"],"non_CDS",rec["ref"],rec["start"],rec["end"],rec["strand"])); used.add(key)
    if not features: raise SystemExit("ERROR: no appropriate CDS or non-CDS RNA features were found in the GFF.")
    if not any(f.id for f in features): raise SystemExit("ERROR: GFF feature identifiers could not be resolved.")
    return features

def subtype(t):
    low=t.lower()
    for canonical in ("rrna","trna","tmrna","ncrna","srna","misc_rna","snrna"):
        if low==canonical: return {"rrna":"rRNA","trna":"tRNA","tmrna":"tmRNA","ncrna":"ncRNA","srna":"sRNA","misc_rna":"misc_RNA","snrna":"snRNA"}[canonical]
    return "other"

def write_csv(path, header, rows):
    with open(path,"w",newline="") as h:
        w=csv.writer(h); w.writerow(header); w.writerows(rows)

def aggregate(sample, features, counts, bins, minimum, cls=None, by_type=False):
    selected=[f for f in features if (cls is None or f.cls==cls) and sum(counts[f])>=minimum]
    groups=defaultdict(list)
    for f in selected: groups[subtype(f.type) if by_type else f.cls].append(f)
    rows=[]
    for group in sorted(groups):
        fs=groups[group]
        for b in range(bins):
            vals=[counts[f][b]/sum(counts[f]) for f in fs]
            rows.append([sample,group,b+1,100*b/bins,100*(b+1)/bins,
                         sum(vals)/len(vals),statistics.median(vals),len(fs),sum(counts[f][b] for f in fs)])
    return rows

def plot_profiles(prefix, bins, cds_rows, nc_rows, type_rows):
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc: raise SystemExit("ERROR: matplotlib is required to generate metagene plots.") from exc
    def one(rows, title, path):
        fig,ax=plt.subplots(figsize=(8,4.5)); x=[(r[3]+r[4])/2 for r in rows]
        if rows: ax.plot(x,[r[5] for r in rows],label="Mean"); ax.plot(x,[r[6] for r in rows],label="Median",alpha=.75)
        ax.set(xlabel="Normalized feature position (0% = 5′ end; 100% = 3′ end)",ylabel="Normalized read fraction",title=title,xlim=(0,100)); ax.legend(); fig.tight_layout(); fig.savefig(path,dpi=150); plt.close(fig)
    one(cds_rows,"CDS metagene profile",prefix+".metagene_CDS_profile.png")
    one(nc_rows,"Non-CDS RNA metagene profile",prefix+".metagene_non_CDS_profile.png")
    fig,ax=plt.subplots(figsize=(8,4.5))
    for rows,label in ((cds_rows,"CDS"),(nc_rows,"non-CDS")):
        if rows: ax.plot([(r[3]+r[4])/2 for r in rows],[r[5] for r in rows],label=label)
    ax.set(xlabel="Normalized feature position (0% = 5′ end; 100% = 3′ end)",ylabel="Mean normalized read fraction",title="CDS vs non-CDS",xlim=(0,100)); ax.legend(); fig.tight_layout(); fig.savefig(prefix+".metagene_CDS_vs_non_CDS.png",dpi=150); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for typ in sorted({r[1] for r in type_rows}):
        rr=[r for r in type_rows if r[1]==typ]; ax.plot([(r[3]+r[4])/2 for r in rr],[r[5] for r in rr],label=typ)
    ax.set(xlabel="Normalized feature position (0% = 5′ end; 100% = 3′ end)",ylabel="Mean normalized read fraction",title="Non-CDS profiles by type",xlim=(0,100));
    if type_rows: ax.legend(); fig.tight_layout(); fig.savefig(prefix+".metagene_non_CDS_by_type.png",dpi=150); plt.close(fig)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--bam",required=True); p.add_argument("--gff",required=True); p.add_argument("--gene-id",required=True); p.add_argument("--bins",type=int,default=100); p.add_argument("--min-feature-reads",type=int,default=10); p.add_argument("--output-prefix",required=True); p.add_argument("--paired",action="store_true"); a=p.parse_args()
    for path,label in ((a.bam,"BAM"),(a.gff,"GFF")):
        if not os.path.isfile(path): p.error(f"{label} is missing: {path}")
    if not (os.path.isfile(a.bam+".bai") or os.path.isfile(os.path.splitext(a.bam)[0]+".bai")): p.error(f"BAM index is missing for {a.bam}")
    if a.bins<1 or a.min_feature_reads<1: p.error("--bins and --min-feature-reads must be positive")
    features=parse_gff(a.gff,a.gene_id); sample=os.path.basename(a.output_prefix)
    bam=pysam.AlignmentFile(a.bam,"rb"); common=set(bam.references)&{f.ref for f in features}
    if not common: raise SystemExit("ERROR: BAM and GFF reference names do not match.")
    features=[f for f in features if f.ref in common]; by_ref=defaultdict(list)
    for f in features: by_ref[f.ref].append(f)
    counts={f:[0]*a.bins for f in features}; unique=Counter(); ambiguous=Counter(); type_unique=Counter(); total=outside=0; detail=[]
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
        if a.paired:
            if not read.is_paired or not read.is_read1 or read.mate_is_unmapped or read.next_reference_id!=read.reference_id or read.template_length==0: continue
            left=min(read.reference_start,read.next_reference_start); right=left+abs(read.template_length); midpoint=(left+right)/2; ident=read.query_name
        else:
            if read.reference_end is None: continue
            midpoint=(read.reference_start+read.reference_end)/2; ident=read.query_name
        total+=1; hits=[f for f in by_ref[read.reference_name] if f.start<=midpoint<f.end]
        any_hit=False
        for cls in ("CDS","non_CDS"):
            match=[f for f in hits if f.cls==cls]
            if len(match)==1:
                any_hit=True; f=match[0]; pos=100*((midpoint-f.start)/(f.end-f.start) if f.strand=="+" else (f.end-midpoint)/(f.end-f.start)); pos=max(0,min(100,pos)); b=min(a.bins-1,int(pos*a.bins/100)); counts[f][b]+=1; unique[cls]+=1
                if cls=="non_CDS": type_unique[subtype(f.type)]+=1
                detail.append([sample,f.id,f.name,f.type,f.cls,f.ref,f.start+1,f.end,f.strand,ident,midpoint,pos,read.mapping_quality,"unique"])
            elif len(match)>1:
                any_hit=True; ambiguous[cls]+=1
                for f in match: detail.append([sample,f.id,f.name,f.type,f.cls,f.ref,f.start+1,f.end,f.strand,ident,midpoint,"",read.mapping_quality,"ambiguous_excluded"])
        if not any_hit: outside+=1
    bam.close()
    write_csv(a.output_prefix+".feature_read_positions.csv",["sample","feature_id","feature_name","feature_type","feature_class","reference","feature_start","feature_end","strand","read_or_fragment_id","alignment_midpoint","position_percent","mapping_quality","assignment_status"],detail)
    binrows=[]
    for f in features:
        n=sum(counts[f])
        for b,c in enumerate(counts[f]): binrows.append([sample,f.id,f.name,f.type,f.cls,b+1,100*b/a.bins,100*(b+1)/a.bins,c,c/n if n else 0])
    write_csv(a.output_prefix+".feature_position_bins.csv",["sample","feature_id","feature_name","feature_type","feature_class","bin","bin_start_percent","bin_end_percent","read_count","fraction_of_feature_reads"],binrows)
    cds=aggregate(sample,features,counts,a.bins,a.min_feature_reads,"CDS"); nc=aggregate(sample,features,counts,a.bins,a.min_feature_reads,"non_CDS"); types=aggregate(sample,[f for f in features if f.cls=="non_CDS"],counts,a.bins,a.min_feature_reads,by_type=True)
    header=["sample","feature_class","bin","bin_start_percent","bin_end_percent","mean_fraction","median_fraction","features_contributing","raw_read_count"]
    write_csv(a.output_prefix+".metagene_CDS_profile.csv",header,cds); write_csv(a.output_prefix+".metagene_non_CDS_profile.csv",header,nc)
    write_csv(a.output_prefix+".metagene_non_CDS_by_type.csv",["sample","feature_type"]+header[2:],types)
    summary=[["total_mapped_reads_or_fragments_examined",total],["reads_or_fragments_outside_all_analyzed_features",outside]]
    for cls in ("CDS","non_CDS"):
        fs=[f for f in features if f.cls==cls]; summary += [[f"{cls}.uniquely_assigned",unique[cls]],[f"{cls}.ambiguous",ambiguous[cls]],[f"{cls}.features_with_reads",sum(sum(counts[f])>0 for f in fs)],[f"{cls}.features_in_metagene",sum(sum(counts[f])>=a.min_feature_reads for f in fs)]]
    for typ in sorted({subtype(f.type) for f in features if f.cls == "non_CDS"}):
        fs = [f for f in features if f.cls == "non_CDS" and subtype(f.type) == typ]
        summary.extend([[f"non_CDS.{typ}.annotated_features", len(fs)],
                        [f"non_CDS.{typ}.uniquely_assigned", type_unique[typ]],
                        [f"non_CDS.{typ}.features_with_reads", sum(sum(counts[f]) > 0 for f in fs)],
                        [f"non_CDS.{typ}.features_in_metagene", sum(sum(counts[f]) >= a.min_feature_reads for f in fs)]])
    write_csv(a.output_prefix+".feature_position_summary.csv",["metric","value"],summary); plot_profiles(a.output_prefix,a.bins,cds,nc,types)
if __name__=="__main__": main()
