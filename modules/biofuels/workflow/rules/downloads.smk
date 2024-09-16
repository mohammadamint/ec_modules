# We recommend adding rules that download necessary files here.
if config["use_default_user_resources"]:
    rule download_units:
        message: "Download spatial units."
        params:
            url = internal["resources"]["default_user_shapes"],
        output: "resources/user/spatial_units.geojson"
        conda: "../envs/shell.yaml"
        localrule: True
        shell: "curl -sSLo {output} '{params.url}'"

rule download_potentials:
    message: "Download potential data."
    params: url = internal["resources"]["potentials"]
    output: protected("resources/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"

rule unzip_potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    output:
        population = "resources/{resolution}/population.csv",
        land_cover = "resources/{resolution}/land-cover.csv"
    conda: "../envs/shell.yaml"
    shell: "unzip {input} '{wildcards.resolution}/*' -d resources"

rule download_biofuel_potentials_and_costs:
    message: "Download raw biofuel potential and cost data."
    params: url = internal["resources"]["biofuel-potentials-and-costs"]
    output: protected("data/automatic/raw-biofuel-potentials-and-costs.xlsx")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"