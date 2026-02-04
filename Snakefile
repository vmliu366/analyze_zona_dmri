from snakebids import bids
from snakebids import generate_inputs
from pathlib import Path

configfile: 'config.yml'


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
        inputs['dwi'].expand(
        bids(suffix='dwi.mif',datatype='dwi',**inputs['dwi'].wildcards))

rule import_mif:
    input:
        dwi=inputs['dwi'].path,
        bval=sidecar(inputs['dwi'].path,'.bval'),
        bvec=sidecar(inputs['dwi'].path,'.bvec'),
    output:
        dwi=bids(root='results',suffix='dwi.mif',datatype='dwi',**inputs['dwi'].wildcards)
    shell:
        'mrconvert {input.dwi} -fslgrad {input.bvec} {input.bval} {output.dwi}'



