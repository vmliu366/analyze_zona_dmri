from snakebids import bids
from snakebids import generate_inputs
from pathlib import Path

configfile: 'config.yml'

root='results'

def sidecar(path: Path, new_suffix: str) -> Path:
    """
    Replace .nii or .nii.gz with a new suffix (e.g. .bval, .bvec, .json).

    Examples
    --------
    sidecar(Path("dwi.nii.gz"), ".bval") -> Path("dwi.bval")
    sidecar(Path("dwi.nii"), ".json")    -> Path("dwi.json")
    """
    path = Path(path)

    if path.suffixes[-2:] == [".nii", ".gz"]:
        stem = path.name[:-7]   # remove ".nii.gz"
    else:
        stem = path.stem        # removes single suffix

    return path.with_name(stem + new_suffix)



inputs = generate_inputs(bids_dir=config['bids_dir'], pybids_inputs=config['pybids_inputs'])

rule all:
    input:
        metrics =inputs['dwi'].expand(
        bids(root=root,suffix='{metric}.nii',datatype='dwi',**inputs['dwi'].wildcards),
            metric=['FA','ADC','AD','RD','V1','V2','V3']),
        tracks=inputs['dwi'].expand(
        bids(root=root,suffix='tracks.tck',algo='{algo}',select='{select}',datatype='dwi',**inputs['dwi'].wildcards),
            select=config['tracking']['select'],
            algo=config['tracking']['algo'])

rule import_mif:
    input:
        dwi=inputs['dwi'].path,
        bval=sidecar(inputs['dwi'].path,'.bval'),
        bvec=sidecar(inputs['dwi'].path,'.bvec'),
    output:
        dwi=bids(root=root,suffix='dwi.mif',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'mrconvert {input.dwi} -fslgrad {input.bvec} {input.bval} {output.dwi}'

rule dwi2tensor:
    input:
        dwi=bids(root=root,suffix='dwi.mif',datatype='dwi',**inputs['dwi'].wildcards)
    output:
        dt=bids(root=root,suffix='dt.mif',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'dwi2tensor {input} {output}'


rule tensor2metric_scalars:
    input:
        dt=bids(root=root,suffix='dt.mif',datatype='dwi',**inputs['dwi'].wildcards)
    params:
        opts=lambda wildcards: '-{metric_lower}'.format(metric_lower=str(wildcards.metric).lower())
    output:
        metric=bids(root=root,suffix='{metric,FA|ADC|AD|RD}.nii',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'tensor2metric {input} {params} {output}'


rule tensor2eigenvec:
    input:
        dt=bids(root=root,suffix='dt.mif',datatype='dwi',**inputs['dwi'].wildcards)
    output:
        metric=bids(root=root,suffix='V{num}.nii',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'tensor2metric {input} -num {wildcards.num} -modulate none -vector {output}'

rule get_fa_mask:
    input:
        metric=bids(root=root,suffix='FA.nii',datatype='dwi',**inputs['dwi'].wildcards)
    output:
        mask=bids(root=root,suffix='mask.nii',desc='thFA0p{th}',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'c3d {input} -threshold 0.{wildcards.th} 1 1 0 -o {output.mask}'
    

rule tracking:
    input:
        dwi=bids(root=root,suffix='dwi.mif',datatype='dwi',**inputs['dwi'].wildcards),
        seed=bids(root=root,suffix='mask.nii',desc='thFA0p5',datatype='dwi',**inputs['dwi'].wildcards)
    params:
        algo=lambda wildcards: 'Tensor_Det' if wildcards.algo == 'det' else 'Tensor_Prob'
    output:
        tck=bids(root=root,suffix='tracks.tck',algo='{algo,prob|det}',select='{select}',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'tckgen {input.dwi} -algorithm {params.algo} -select {wildcards.select} '
        '-seed_image {input.seed} '
        '{output.tck}'

