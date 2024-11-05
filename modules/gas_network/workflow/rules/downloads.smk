"""Rules to download necessary resources."""

if config["use_default_user_resources"]:
    rule user_input_shapes:
        message:
            "Download spatial zones."
        params:
            url=internal["resources"]["spatial_units"],
        output:
            "resources/user/shapes_e_highways.geojson",
        # conda:
        #     "../envs/shell.yaml"
        localrule: True
        shell:
            "curl -sSLo {output} '{params.url}'"   

rule download_sci_grid_data:
    message: "Download gas infrastructure data from SciGRID_gas IGGIELGN"
    params:
        url = internal["resources"]["SciGRID_gas"]
    output: "resources/automatic/gas_grid.zip"
    # conda: "../envs/shell.yaml"
    shell:
        """curl -sSLo {output} {params.url}"""

