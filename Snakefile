from snakebids import bids
from snakebids import generate_inputs
from pathlib import Path


configfile: "config.yml"

inputs = generate_inputs(
    bids_dir=config["bids_dir"], pybids_inputs=config["pybids_inputs"],
    derivatives=config["derivatives"],
)



root = "results"

NODE_START = list(map(int, config["connectome"]["node_start"]))
NODE_END   = list(map(int, config["connectome"]["node_end"]))

EDGE_PAIRS = [(s, e) for s in NODE_START for e in NODE_END]

EDGEPAIR_STRS = [f"{a}-{b}" for a, b in EDGE_PAIRS]
EDGEPAIR_DIR_TMPL = bids(
    root=root,
    suffix="edgepairs",
    res="{res}",
    seedmask="{seedmask}",
    algo="{algo}",
    select="{select}",
    nodes="{nodes}",
    datatype="dwi",
    **inputs["dwi"].wildcards,
)
EDGEPAIR_TCK_TMPL = str(Path(EDGEPAIR_DIR_TMPL) / "edge-{edgepair}.tck")


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



rule all:
    input:
        metrics=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="{metric}.nii",
                datatype="dwi",
                res="{res}",
                **inputs["dwi"].wildcards,
            ),
            metric=["FA", "ADC", "AD", "RD", "V1", "V2", "V3"],
            res=config["downsample_res"],
        ),
        tracks=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="tracks.tck",
                res="{res}",
                seedmask="{seedmask}",
                algo="{algo}",
                select="{select}",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            res=config["downsample_res"],
            seedmask=config["tracking"]["seedmask"],
            select=config["tracking"]["select"],
            algo=config["tracking"]["algo"],
        ),
        exemplars=inputs["dwi"].expand(
            bids(
            root=root,
            suffix="exemplars.tck",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            nodes="{nodes}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
            res=config["downsample_res"],
            seedmask=config["tracking"]["seedmask"],
            select=config["tracking"]["select"],
            algo=config["tracking"]["algo"],
            nodes=config['connectome']['nodes'],

        ),
        edge_tcks=inputs["dwi"].expand(
            EDGEPAIR_TCK_TMPL,
            res=config["downsample_res"],
            seedmask=config["tracking"]["seedmask"],
            select=config["tracking"]["select"],
            algo=config["tracking"]["algo"],
            nodes=config["connectome"]["nodes"],
            edgepair=EDGEPAIR_STRS,
        ),
        meshes=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="dseg.obj",
                datatype="dwi",
                res="{res}",
                **inputs["dwi"].wildcards,
            ),
            nodes=config["connectome"]["nodes"],
            res=config["downsample_res"],
        ),
        pca_labels=inputs["dwi"].expand(
            bids(
                root=root,
                suffix="streamlabels.txt",
                res="{res}",
                seedmask="{seedmask}",
                algo="{algo}",
                select="{select}",
                points="{points}",
                desc="pcakmeans",
                scalargroup="{scalargroup}",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            res=config["downsample_res"],
            seedmask=config["tracking"]["seedmask"],
            select=config["tracking"]["select"],
            algo=config["tracking"]["algo"],
            points=config["clustering"]["points"],
            scalargroup=config["clustering"]["scalargroups"].keys(),
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


rule downsample_dwi:
    input:
        dwi=bids(root=root, suffix="dwi.mif", datatype="dwi", **inputs["dwi"].wildcards),
    params:
        voxel_size=lambda wildcards: str(wildcards.res).replace("p", "."),
    output:
        dwi=bids(
            root=root,
            suffix="dwi.mif",
            res="{res}mm",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "mrgrid {input.dwi} regrid {output.dwi} -voxel {params.voxel_size} -interp cubic"


rule dwi2tensor_unmasked:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.mif",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
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

rule anat2dwi_reg:
    """
    Create transforms from subject anat space (T1w) to subject DWI b0 space
    using antsRegistration (affine-only)
    """
    input:
        fixed=bids(root=root, suffix="b0.nii", datatype="dwi", **inputs["dwi"].wildcards),
        moving=inputs["t1w"].path,
    params:
        out_prefix=bids(
            root=root,
            datatype="dwi",
            desc="from-T1w_to-dwi",
            suffix="",
            extension="",
            **inputs["dwi"].wildcards,
        ),
    output:
        affine=bids(
            root=root,
            datatype="dwi",
            desc="from-T1w_to-dwi",
            suffix="0GenericAffine",
            extension=".mat",
            **inputs["dwi"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=60,
    shell:
        r"""
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}

        mkdir -p "$(dirname "{params.out_prefix}")"
        echo "created $(dirname "{params.out_prefix}")"

        antsRegistrationSyN.sh \
            -d 3 \
            -f {input.fixed} \
            -m {input.moving} \
            -o {params.out_prefix} \
            -t r \
            -n {threads}
        """


rule anat2dwi_apply:
    """
    Apply calculated affine to T1w (Linear) and dseg (MultiLabel)
    """
    input:
        t1w=inputs["t1w"].path,
        dseg=inputs["dseg"].path,
        ref=bids(root=root, suffix="b0.nii", datatype="dwi", **inputs["dwi"].wildcards),
        affine=rules.anat2dwi_reg.output.affine,
    params:
        out_dir=lambda wildcards: str(Path(bids(
            root=root, datatype="anat", space="dwi",
            suffix="T1w", extension=".nii.gz",
            **inputs["dwi"].wildcards
        )).parent),
    output:
        t1w_in_dwi=bids(
            root=root,
            datatype="anat",
            space="dwi",
            suffix="T1w",
            extension=".nii.gz",
            **inputs["dwi"].wildcards,
        ),
        dseg_in_dwi=bids(
            root=root,
            datatype="anat",
            space="dwi",
            suffix="dseg",
            extension=".nii.gz",
            **inputs["dwi"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        time=30,
    shell:
        r"""
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}
        mkdir -p {params.out_dir}
        echo "created {params.out_dir}"

        antsApplyTransforms -v -d 3 -n Linear \
          -i {input.t1w} -r {input.ref} \
          -t {input.affine} \
          -o {output.t1w_in_dwi}

        antsApplyTransforms -v -d 3 -n MultiLabel \
          -i {input.dseg} -r {input.ref} \
          -t {input.affine} \
          -o {output.dseg_in_dwi}
        """

rule brain_mask:
    input:
        t1w_in_dwi=bids(
            root=root,
            datatype="anat",
            space="dwi",
            suffix="T1w",
            extension=".nii.gz",
            **inputs["dwi"].wildcards,
        ),
    output:
        mask=bids(
            root=root,
            suffix="mask.nii",
            desc="brain",
            datatype="dwi",
            # **inputs["dwi"].wildcards,
            **(lambda d: (d.pop("desc", None), d)[1])(dict(inputs["dwi"].wildcards)),
        ),
    shell:
        "c3d {input.t1w_in_dwi} -otsu 3 -binarize -o {output.mask}"


rule downsample_mask:
    input:
        mask=bids(
            root=root,
            suffix="mask.nii",
            desc="brain",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    params:
        voxel_size=lambda wildcards: str(wildcards.res).replace("p", "."),
    output:
        mask=bids(
            root=root,
            suffix="mask.nii",
            res="{res}mm",
            desc="brain",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "mrgrid {input.mask} regrid {output.mask} -voxel {params.voxel_size} -interp nearest"



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

rule dwi2tensor_ds:
    input:
        dwi=bids(root=root, suffix="dwi.mif", res="{res}", datatype="dwi", **inputs["dwi"].wildcards),
        mask=bids(
            root=root,
            suffix="mask.nii",res="{res}", 
            desc="brain",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    output:
        dt=bids(root=root, suffix="dt.mif", res="{res}", datatype="dwi", **inputs["dwi"].wildcards),
    shell:
        "dwi2tensor {input.dwi} -mask {input.mask} {output.dt}"



rule tensor2metric_scalars:
    input:
        dt=bids(root=root, suffix="dt.mif", res="{res}",datatype="dwi", **inputs["dwi"].wildcards),
    params:
        opts=lambda wildcards: "-{metric_lower}".format(
            metric_lower=str(wildcards.metric).lower()
        ),
    output:
        metric=bids(
            root=root,
            res="{res}",
            suffix="{metric,FA|ADC|AD|RD}.nii",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "tensor2metric {input} {params} {output}"


rule tensor2eigenvec:
    input:
        dt=bids(root=root, suffix="dt.mif",
                res="{res}",
                datatype="dwi", **inputs["dwi"].wildcards),
    output:
        metric=bids(
            root=root, suffix="V{num}.nii",
            res="{res}",
            datatype="dwi", **inputs["dwi"].wildcards
        ),
    shell:
        "tensor2metric {input} -num {wildcards.num} -modulate none -vector {output}"


rule get_fa_mask:
    input:
        metric=bids(
            root=root, suffix="FA.nii",res="{res}", 
            datatype="dwi", **inputs["dwi"].wildcards
        ),
    output:
        mask=bids(
            root=root,
            suffix="mask.nii",
            res="{res}", 
            desc="thFA0p{th}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "c3d {input} -threshold 0.{wildcards.th} 1 1 0 -o {output.mask}"


rule tracking:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.mif",
            datatype="dwi",
            res="{res}",
            **inputs["dwi"].wildcards,
        ),
        seed=bids(
            root=root,
            suffix="mask.nii",res="{res}", 
            desc="{seedmask}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        mask=bids(
            root=root,
            suffix="mask.nii",res="{res}", 
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
            res="{res}",
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
            res="{res}",
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
            res="{res}",
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
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        img=bids(
            root=root, suffix="{scalar}.nii", res="{res}",
            datatype="dwi", **inputs["dwi"].wildcards
        ),
    output:
        txt=bids(
            root=root,
            suffix="trackscalars.txt",
            res="{res}",
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
        fa=bids(root=root, suffix="FA.nii", res="{res}",datatype="dwi", **inputs["dwi"].wildcards),
    params:
        coord_pat=bids(
            root=root, suffix="coord%d.nii",res="{res}", datatype="dwi", **inputs["dwi"].wildcards
        ),
    output:
        coords=expand(
            bids(
                root=root,
                suffix="coord{i}.nii",
                res="{res}",
                datatype="dwi",
                **inputs["dwi"].wildcards,
            ),
            allow_missing=True,
            i=[0, 1, 2],
        ),
    shell:
        "c3d {input} -cmv -oo {params.coord_pat}"


rule concat_scalars:
    """produces a N_scalars * N_streamlines * N_points array"""
    input:
        scalar_txts=lambda wildcards: expand(
            bids(
                root=root,
                suffix="trackscalars.txt",
                res="{res}",
                seedmask="{seedmask}",
                algo="{algo,prob|det}",
                select="{select}",
                points="{points}",
                datatype="dwi",
                scalar="{scalar}",
                **inputs["dwi"].wildcards,
            ),
            allow_missing=True,
            scalar=config["clustering"]["scalargroups"][wildcards.scalargroup],
        ),
    output:
        npy=bids(
            root=root,
            suffix="concatscalars.npy",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            scalargroup="{scalargroup}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    run:
        import numpy as np

        stacked = np.stack(
            [np.loadtxt(scalar_txt, skiprows=1) for scalar_txt in input.scalar_txts]
        )
        np.save(output.npy, stacked)


rule pca_kmeans:
    input:
        npy=bids(
            root=root,
            suffix="concatscalars.npy",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo,prob|det}",
            select="{select}",
            points="{points}",
            scalargroup="{scalargroup}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    params:
        n_components=20,
        k=5,
    output:
        npy=bids(
            root=root,
            suffix="streamlabels.txt",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            points="{points}",
            desc="pcakmeans",
            scalargroup="{scalargroup}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    script:
        "scripts/pca_kmeans.py"


#-- dseg labels
# rule reslice_dseg_to_dwi:
#     """ uses existing itksnap rigid transform, may need non-linear for best match"""
#     input:
#         dseg = inputs['dseg'].path,
#         xfm = inputs['xfm_to_dwi'].path,
#         ref=bids(
#             root=root, suffix="FA.nii", res="{res}",
#             datatype="dwi", **inputs["dwi"].wildcards
#         ),

#     output:
#         dseg=bids(
#             root=root, suffix="dseg.nii.gz",desc="{desc}", res="{res}",
#             datatype="dwi", **inputs["dwi"].wildcards
#         ),
#     shell:
#         # "antsApplyTransforms -d 3 -n NearestNeighbor -i {input.dseg} -t {input.xfm} -r {input.ref} -o {output.dseg} -u int"


rule tck2connectome:
    input:
        tck=bids(
            root=root,
            suffix="tracks.tck",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        # dseg=bids(
        #     root=root, suffix="dseg.nii.gz",desc="{desc}", res="{res}",
        #     datatype="dwi", **inputs["dwi"].wildcards
        # ),
        dseg=rules.anat2dwi_apply.output.dseg_in_dwi,
        
    output:
        connectome=bids(
            root=root,
            suffix="connectome.csv",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            nodes="{desc}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        assignments=bids(
            root=root,
            suffix="assignments.txt",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            nodes="{desc}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),

    shell:
        "tck2connectome {input.tck} {input.dseg} {output.connectome} -out_assignments {output.assignments} -assignment_radial_search 1"



# rule connectome2tck:
#     input:
#         tck=bids(
#             root=root,
#             suffix="tracks.tck",
#             res="{res}",
#             seedmask="{seedmask}",
#             algo="{algo}",
#             select="{select}",
#             datatype="dwi",
#             **inputs["dwi"].wildcards,
#         ),
#         assignments=bids(
#             root=root,
#             suffix="assignments.txt",
#             res="{res}",
#             seedmask="{seedmask}",
#             algo="{algo}",
#             select="{select}",
#             nodes="{desc}",
#             datatype="dwi",
#             **inputs["dwi"].wildcards,
#         ),
#     params:
#         bundle_prefix = lambda wildcards, output: str(Path(output.bundle_dir) / 'bundles_')
#     output:
#         bundle_dir=directory(bids(
#             root=root,
#             suffix="bundles",
#             res="{res}",
#             seedmask="{seedmask}",
#             algo="{algo}",
#             select="{select}",
#             nodes="{desc}",
#             datatype="dwi",
#             **inputs["dwi"].wildcards,
#         )),
#     shell:
#         "mkdir -p {output.bundle_dir} && "
#         "connectome2tck {input.tck} {input.assignments} {params.bundle_prefix}"

rule connectome2tck_edgepair:
    input:
        tck=bids(
            root=root, suffix="tracks.tck",
            res="{res}", seedmask="{seedmask}", algo="{algo}", select="{select}",
            datatype="dwi", **inputs["dwi"].wildcards,
        ),
        assignments=bids(
            root=root, suffix="assignments.txt",
            res="{res}", seedmask="{seedmask}", algo="{algo}", select="{select}",
            nodes="{desc}",
            datatype="dwi", **inputs["dwi"].wildcards,
        ),
    params:
        nodes=lambda wc: wc.edgepair.replace("-", ","),
        out_dir=bids(
            root=root, suffix="edgepairs",
            res="{res}", seedmask="{seedmask}", algo="{algo}", select="{select}",
            nodes="{desc}",
            datatype="dwi", **inputs["dwi"].wildcards,
        ),
    output:
        tck=str(Path(
            bids(
                root=root, suffix="edgepairs",
                res="{res}", seedmask="{seedmask}", algo="{algo}", select="{select}",
                nodes="{desc}",
                datatype="dwi", **inputs["dwi"].wildcards,
            )
        ) / "edge-{edgepair}.tck"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        connectome2tck {input.tck} {input.assignments} {output.tck} \
          -nodes {params.nodes} -exclusive -files single
        """



rule connectome2tck_exemplars:
    input:
        tck=bids(
            root=root,
            suffix="tracks.tck",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        assignments=bids(
            root=root,
            suffix="assignments.txt",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            nodes="{desc}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
        # dseg=bids(
        #     root=root, suffix="dseg.nii.gz",desc="{desc}", res="{res}",
        #     datatype="dwi", **inputs["dwi"].wildcards
        # ),
        dseg=rules.anat2dwi_apply.output.dseg_in_dwi,

    output:
        exemplars=bids(
            root=root,
            suffix="exemplars.tck",
            res="{res}",
            seedmask="{seedmask}",
            algo="{algo}",
            select="{select}",
            nodes="{desc}",
            datatype="dwi",
            **inputs["dwi"].wildcards,
        ),
    shell:
        "connectome2tck {input.tck} {input.assignments} {output.exemplars} "
        " -files single -exemplars {input.dseg}"
        
rule label2mesh:
    input:
        # dseg=bids(
        #     root=root, suffix="dseg.nii.gz",desc="{desc}", res="{res}",
        #     datatype="dwi", **inputs["dwi"].wildcards
        # ),
        dseg=rules.anat2dwi_apply.output.dseg_in_dwi,
    output:
        dseg=bids(
            root=root, suffix="dseg.obj", res="{res}",
            datatype="dwi", **inputs["dwi"].wildcards
        ),
    shell:
        "label2mesh {input} {output}"


