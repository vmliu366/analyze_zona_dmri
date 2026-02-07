from snakebids import bids
from snakebids import generate_inputs
from pathlib import Path


configfile: "config.yml"


root = "results"


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
        stem = path.name[:-7]  # remove ".nii.gz"
    else:
        stem = path.stem  # removes single suffix

    return path.with_name(stem + new_suffix)


inputs = generate_inputs(
    bids_dir=config["bids_dir"], pybids_inputs=config["pybids_inputs"]
)


rule all:
    input:
        metrics=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="{metric}.nii",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            metric=["FA", "ADC", "AD", "RD", "V1", "V2", "V3"],
        ),
        tracks=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="tracks.tck",
                seedmask="{seedmask}",
                algo="{algo}",
                select="{select}",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            seedmask=config["tracking"]["seedmask"],
            select=config["tracking"]["select"],
            algo=config["tracking"]["algo"],
        ),
        scalars_npy=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="concatscalars.npy",
                seedmask="{seedmask}",
                algo="{algo}",
                select="{select}",
                points="{points}",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            seedmask=config["tracking"]["seedmask"],
            select=config["tracking"]["select"],
            algo=config["tracking"]["algo"],
            points=config["clustering"]["points"],
        ),


rule import_mif:
    input:
        dwi=inputs["dwi"].path,
        bval=sidecar(inputs["dwi"].path, ".bval"),
        bvec=sidecar(inputs["dwi"].path, ".bvec"),
    output:
        dwi=bids(root=root, suffix="dwi.mif", datatype="dwi", **inputs["dwi"].wildcards),
    shell:
        "mrconvert {input.dwi} -fslgrad {input.bvec} {input.bval} {output.dwi}"


rule dwi2tensor_unmasked:
    input:
        dwi=bids(root=root, suffix="dwi.mif", datatype="dwi", **inputs["dwi"].wildcards),
    output:
        b0=bids(root=root, suffix="b0.nii", datatype="dwi", **inputs["dwi"].wildcards),
        dt=temp(
            bids(
                root=root,
                suffix="dtunmasked.mif",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            )
        ),
    shell:
        "dwi2tensor {input} -b0 {output.b0} {output.dt}"


rule brain_mask:
    input:
        b0=bids(root=root, suffix="b0.nii", datatype="dwi", **inputs["dwi"].wildcards),
    output:
        mask=bids(
            root=root,
            suffix="mask.nii",
            desc="brain",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "c3d {input} -otsu 3 -binarize {output}"


rule dwi2tensor:
    input:
        dwi=bids(root=root, suffix="dwi.mif", datatype="dwi", **inputs["dwi"].wildcards),
        mask=bids(
            root=root,
            suffix="mask.nii",
            desc="brain",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    output:
        dt=bids(root=root, suffix="dt.mif", datatype="dwi", **inputs["dwi"].wildcards),
    shell:
        "dwi2tensor {input.dwi} -mask {input.mask} {output.dt}"


rule tensor2metric_scalars:
    input:
        dt=bids(root=root, suffix="dt.mif", datatype="dwi", **inputs["dwi"].wildcards),
    params:
        opts=lambda wildcards: "-{metric_lower}".format(
            metric_lower=str(wildcards.metric).lower()
        ),
    output:
        metric=bids(
            root=root,
            suffix="{metric,FA|ADC|AD|RD}.nii",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "tensor2metric {input} {params} {output}"


rule tensor2eigenvec:
    input:
        dt=bids(root=root, suffix="dt.mif", datatype="dwi", **inputs["dwi"].wildcards),
    output:
        metric=bids(
            root=root, suffix="V{num}.nii", datatype="dwi", **inputs["dwi"].wildcards
        ),
    shell:
        "tensor2metric {input} -num {wildcards.num} -modulate none -vector {output}"


rule get_fa_mask:
    input:
        metric=bids(
            root=root, suffix="FA.nii", datatype="dwi", **inputs["dwi"].wildcards
        ),
    output:
        mask=bids(
            root=root,
            suffix="mask.nii",
            desc="thFA0p{th}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "c3d {input} -threshold 0.{wildcards.th} 1 1 0 -o {output.mask}"


rule tracking:
    input:
        dwi=bids(root=root, suffix="dwi.mif", datatype="dwi", **inputs["dwi"].wildcards),
        seed=bids(
            root=root,
            suffix="mask.nii",
            desc="{seedmask}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        mask=bids(
            root=root,
            suffix="mask.nii",
            desc="brain",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    params:
        algo=lambda wildcards: (
            "Tensor_Det" if wildcards.algo == "det" else "Tensor_Prob"
        ),
    output:
        tck=bids(
            root=root,
            suffix="tracks.tck",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "tckgen {input.dwi} -algorithm {params.algo} -select {wildcards.select} "
        "-seed_image {input.seed} "
        "-mask {input.mask} "
        "{output.tck}"


rule resample_tck:
    input:
        tck=bids(
            root=root,
            suffix="tracks.tck",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    output:
        tck=bids(
            root=root,
            suffix="resampledtracks.tck",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "tckresample {input} -num_points {wildcards.points} {output}"


rule sample_scalars_on_tck:
    input:
        tck=bids(
            root=root,
            suffix="resampledtracks.tck",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        img=bids(
            root=root, suffix="{scalar}.nii", datatype="dwi", **inputs["dwi"].wildcards
        ),
    output:
        txt=bids(
            root=root,
            suffix="trackscalars.txt",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            datatype="dwi",
            scalar="{scalar}",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "tcksample {input.tck} {input.img} {output.txt}"


rule make_coords_vol:
    """ x y z coords for mapping coordinate to tck scalar file. note: it's already in .tck itself, but 
    having them as volumes  then scalars just makes the workflow easier"""
    input:
        b0=bids(root=root, suffix="b0.nii", datatype="dwi", **inputs["dwi"].wildcards),
    params:
        coord_pat=bids(
            root=root, suffix="coord%d.nii", datatype="dwi", **inputs["dwi"].wildcards
        ),
    output:
        coords=expand(
            bids(
                root=root,
                suffix="coord{i}.nii",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            allow_missing=True,
            i=[0, 1, 2],
        ),
    shell:
        "c3d {input.b0} -cmv -oo {params.coord_pat}"


rule concat_scalars:
    input:
        scalar_txts=expand(
            bids(
                root=root,
                suffix="trackscalars.txt",
                seedmask="{seedmask}",
                algo="{algo,prob|det}",
                select="{select}",
                points="{points}",
                datatype="dwi",
                scalar="{scalar}",
                **inputs["dwi"].wildcards,
            ),
            allow_missing=True,
            scalar=config["clustering"]["scalars"],
        ),
    output:
        npy=bids(
            root=root,
            suffix="concatscalars.npy",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    run:
        import numpy as np

        stacked = np.stack(
            [np.loadtxt(scalar_txt, skiprows=1) for scalar_txt in input.scalar_txts]
        )
        print(stacked.shape)
        np.save(output.npy, stacked)
