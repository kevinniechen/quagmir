{
  "id": "https://api.sbgenomics.com/v2/apps/nikola_tesic/quagmir-dev/quagmir/53/raw/",
  "class": "CommandLineTool",
  "label": "QuagmiR",
  "description": "**QuagmiR** is a python tool for mapping of miRNA sequences and detection and quantification of different isomiRs. For mapping, the tool searches for exact match in the conserved middle part of the miRNA sequence, and once the match is found it uses filtering based on other information to decide whether to discard the sequence or not. It uses FASTQ files with miRNA-Seq reads as input.\n\n**QuagmiR** requires the following input files and parameters:\n\n* **Input reads**: MiRNA-Seq reads in FASTQ format (gzipped input accepted). If multiple files are provided, all available threads on the instance will be used so the task's run time is made as short as possible by running **QuagmiR** on samples in parallel.\n\n* **Motif-consensus file**: Reference file containing names of miRNA sequences, their conserved parts and miRNA sequences themselves. For this, [motif_list_hsa.fa](https://cgc.sbgenomics.com/public/files/5a8dab554f0c3a81d19fca69/) can be used, which contains human miRNAs. Otherwise, any list of miRNAs can be used as **Reference file** as long as it adheres to the format of [motif_list_hsa.fa](https://cgc.sbgenomics.com/public/files/5a8dab554f0c3a81d19fca69/) (for example, if there is interest only in specific miRNAs, only those can be used as reference, and if FASTQ files belong to some non-human species, then reference file for that species is needed).\n\n* **Reference file**: Reference TSV containing additional information about miRNA references. Example file can be found on QuagmiR's github repo.\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the end of the page.*\n\n### Common Use Cases\n\nMicroRNAs (abbreviated miRNAs) are small non-coding RNA molecules (containing about 22 nucleotides) found in plants, animals and some viruses, that have a function in RNA silencing and post-transcriptional regulation of gene expression. IsomiRs are sequences that have variations with respect to the reference MiRNA sequence.\n\nEach miRNA is considered to have a conserved part in the middle that is the same for all of its isomiRs. **QuagmiR** searches for that part in the input FASTQ file, and than chooses reads among those that pass the filtering thresholds (such as 3' and 5' ends filtering where the ends need to have smaller weighted Levenstein distance from the miRNA than the user defined threshold (which can be switched off altogether)).\n\nAfter filtering, for sequences that are aligned to multiple miRNA reads, the best alignment is chosen based on the total weighted Levenstein distance between them. This can be changed by setting the **Destructive motif pull** parameter to False. That way, all possible alignments will be outputted.\n\nAdditionally, as previously mentioned, different **Reference file**s can be chosen, so that they better fit a specific use-case (i.e. a subset of miRNA references contained in the **Reference file** we provide on the platform, so that only particular types of miRNA that are of interest are explored).\n\nAll these parameters, and several others that are described at the end of this page, provide the user with a wide range of possibilities to use this tool in. The tool can be used for general miRNA mapping and isomiR quantification, but it can be used for more specific needs as well by applying different filters, using different sets of miRNA references, turning the **Destructive motif pull** option on or off, etc.\n\n### Changes Introduced by Seven Bridges\n\n* Collapsed files can be used as an optional input, so for each FASTQ file that has it's COLLAPSED file, **QuagmiR** will skip the collapsing step, so the execution time can be shortened. This is useful when analyzing same input files with different tool options such as filtering thresholds, so the files are collapsed only the first time. If some of the input files have COLLAPSED files and some don't, COLLAPSED files will be created and outputed for the new input files.\n\n### Common Issues and Important Notes\n\n* At the moment, only **Reference file**s for human miRNA have been built and tested, so the tool can't be used on other species without first creating the **Reference file**, which requires information regarding that species' miRNA that may not be known yet (most notably, length and positions of conserved parts within the known miRNA sequences).\n\n### Performance Benchmarking\n\nIn the following table you can find estimates of running times and costs. All input files are FASTQ files containing miRNA-Seq reads. If not specified otherwise, number of threads for the tool will be set to be the same as number of input files, or 36 in case the number of files is greater than 36. The instance is set to either default c4.2xlarge, or, in case more threads are required, to the cheapest instance that has sufficient number of threads.\n\nExecution time depends on many factors, but the most important ones are size of input files and whether the **Ambiguous letters support** is turned on or not.\n\n*Cost can be significantly reduced by **spot instance** usage. Visit [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*         \n\n#####Example runs:\n\n| Experiment type | Input size | Ambiguous letters | Read length | Duration | Cost | Instance (AWS)| Instance type |\n|-----------------------|------------|-----------------|------------|-----------------|-------------|--------------|------------------|-------------|---------------|------------|\n| miRNA-Seq         | 1    | 864 MB     | No        | 15-30         | 17m   | $0.16            | c4.2xlarge      | On-demand        |\n| miRNA-Seq         | 1    | 864 MB     | Yes        | 15-30         | 1h 19m   | $0.73            | c4.2xlarge      | On-demand        |\n| miRNA-Seq         | 72     | 599MB - 872MB      | Yes        | 15-30         | 2h 19m   | $1.24            | c4.8xlarge      | On-demand        |\n\n### API Python Implementation\n\nThe tool's draft task can also be submitted via the **API**. In order to learn how to get your **Authentication token** and **API endpoint** for corresponding platform visit our [documentation](https://github.com/sbg/sevenbridges-python#authentication-and-configuration).\n\n```python\nfrom sevenbridges import Api\n\nauthentication_token, api_endpoint = \"enter_your_token\", \"enter_api_endpoint\"\napi = Api(token=authentication_token, url=api_endpoint)\n# Get project_id/workflow_id from your address bar. Example: https://igor.sbgenomics.com/u/your_username/project/tool\nproject_id, tool_id = \"your_username/project\", \"your_username/project/tool\"\n# Get file names from files in your project. Example: Names are taken from Data/Public Reference Files.\ninputs = {\n    'fastq_files': list(api.files.query(project=project_id, names=['miRNA_test_file.fastq'])), \n    'motif_consensus_file': api.files.query(project=project_id, names=['motif_list_hsa.fa'])[0]\n    'reference_file': api.files.query(project=project_id, names=['rmiRBase21-master.tsv'])[0]\n}\ntask = api.tasks.create(name='QuagmiR - API Example', project=project_id, app=tool_id, inputs=inputs, run=True)\n```\n\nInstructions for installing and configuring the API Python client, are provided on [github](https://github.com/sbg/sevenbridges-python#installation). For more information about using the API Python client, consult [sevenbridges-python documentation](http://sevenbridges-python.readthedocs.io/en/latest/). **More examples** are available [here](https://github.com/sbg/okAPI).\n\nAdditionally, [API R](https://github.com/sbg/sevenbridges-r) and [API Java](https://github.com/sbg/sevenbridges-java) clients are available. To learn more about using these API clients please refer to the [API R client documentation](https://sbg.github.io/sevenbridges-r/), and [API Java client documentation](https://docs.sevenbridges.com/docs/java-library-quickstart).\n\n### References\n\n[1] [Github repo](https://github.com/duxan/quagmir)\n\n[2] Authors and their contact: [Xavier Bofill De Ross](https://github.com/xbdr86), [Kevin Chen](https://github.com/kevchn), [Nikola Tesic](https://github.com/nikola-tesic), [Dusan Randjelovic](https://github.com/duxan)",
  "requirements": [
    {
      "fileDef": [
        {
          "fileContent": {
            "script": "{\n  var motif_consensus_file = $job.inputs.motif_consensus_file.path\n  var reference_file = $job.inputs.reference_file.path\n  var min_ratio = $job.inputs.min_ratio || -1\n  var min_read = $job.inputs.min_read || -1\n  var destructive_motif_pull = $job.inputs.destructive_motif_pull || \"False\"\n  if(destructive_motif_pull == true)destructive_motif_pull = \"True\"\n  var ambiguous_letters = $job.inputs.ambiguous_letters || \"False\"\n  if(ambiguous_letters == true)ambiguous_letters = \"True\"\n  var display_summary = $job.inputs.display_summary || \"False\"\n  if($job.inputs.display_summary == false) var display_summary = \"False\"\n  else var display_summary = \"True\"\n  if($job.inputs.display_sequence_info == false) var display_sequence_info = \"False\"\n  else var display_sequence_info = \"True\"  \n  if($job.inputs.display_nucleotide_dist == false) var display_nucleotide_dist = \"False\"\n  else var display_nucleotide_dist = \"True\"\n  if($job.inputs.display_group_output == false) var display_group_output = \"False\"\n  else var display_group_output = \"True\"\n  if($job.inputs.group_output_name == null ||\n     $job.inputs.group_output_name == \"\") var group_output_name = \"group_output\"\n  else var group_output_name = $job.inputs.group_output_name\n  var mirbase_version = $job.inputs.mirbase_version || \"21\"\n  var source_ontology = $job.inputs.source_ontology || \"miRBase v21 doi:10.25504/fairsharing.hmgte8\"\n  \n  var edit_distance_3p = $job.inputs.edit_distance_treshold_3p || -1\n  var edit_distance_5p = $job.inputs.edit_distance_treshold_5p || -1\n  \n  //edit distance parameters\n  var deletion_score = $job.inputs.deletion_score || 1\n  var insertion_score = $job.inputs.insertion_score || 1\n  if($job.inputs.substitution_score){\n    var ac = $job.inputs.substitution_score.ac || 1\n    var ag = $job.inputs.substitution_score.ag || 1\n    var at = $job.inputs.substitution_score.at || 1\n    var ca = $job.inputs.substitution_score.ca || 1\n    var cg = $job.inputs.substitution_score.cg || 1\n    var ct = $job.inputs.substitution_score.ct || 1\n    var ga = $job.inputs.substitution_score.ga || 1\n    var gc = $job.inputs.substitution_score.gc || 1\n    var gt = $job.inputs.substitution_score.gt || 1\n    var ta = $job.inputs.substitution_score.ta || 1\n    var tc = $job.inputs.substitution_score.tc || 1\n    var tg = $job.inputs.substitution_score.tg || 1\n  }\n  else{\n    var ac = 1\n    var ag = 1\n    var at = 1\n    var ca = 1\n    var cg = 1\n    var ct = 1\n    var ga = 1\n    var gc = 1\n    var gt = 1\n    var ta = 1\n    var tc = 1\n    var tg = 1\n  }\n  \n  return(\"# FILTERING\\n\" + \n  \"min_ratio: \" + min_ratio + \"\\n\" + //.1 \n  \"min_read: \" + min_read + \"\\n\" + //9\n  \"edit_distance_3p: \" + edit_distance_3p + \"\\n\" + //-1\n  \"edit_distance_5p: \" + edit_distance_5p + \"\\n\\n\" + //-1\n\n  \"# DISPLAY\\n\" +\n  \"display_summary: \" + display_summary + \"\\n\" +//True\n  \"display_sequence_info: \" + display_sequence_info + \"\\n\" + //True\n  \"display_nucleotide_dist: \" + display_nucleotide_dist + \"\\n\" + //True\n  \"display_group_output: \" + display_group_output + \"\\n\\n\" + //True\n         \n  \"# AMBIGUOUS LETTERS\\n\" +\n  \"ambiguous_letters: \" + ambiguous_letters + \"\\n\\n\" + //True\n         \n  \"# DESTRUCTIVE MOTIF PULL\\n\" +\n  \"destructive_motif_pull: \" + destructive_motif_pull + \"\\n\\n\" + //False\n\n  \"# INPUT\\n\" +\n  \"data: data\\n\" + // input folder\n  \"motif_consensus_file: \" + motif_consensus_file + \"\\n\" +\n  \"reference_file: \" + reference_file + '\\n' +\n  \"mirbase_version: \" + mirbase_version + \"\\n\" + // version of miRBase being used for reference files\n  \"source_ontology: \" + source_ontology + \"\\n\\n\" +\n  \n  \"# OUTPUT\\n\" +\n  \"group_output_name: \" + group_output_name + \"\\n\\n\" +\n  \n  \"# DISTANCE METRIC\\n\" +\n  \"deletion_score: \" + deletion_score + \"\\n\" +\n  \"insertion_score: \" + insertion_score + \"\\n\" +\n  \"substitution_AG: \" + ag + \" # transition\\n\" +\n  \"substitution_GA: \" + ga + \" # transition\\n\" +\n  \"substitution_CT: \" + ct + \" # transition\\n\" +\n  \"substitution_TC: \" + tc + \" # transition\\n\" +\n  \"substitution_AT: \" + at + \" # transversion\\n\" +\n  \"substitution_TA: \" + ta + \" # transversion\\n\" +\n  \"substitution_AC: \" + ac + \" # transversion\\n\" +\n  \"substitution_CA: \" + ca + \" # transversion\\n\" +\n  \"substitution_GC: \" + gc + \" # transversion\\n\" +\n  \"substitution_CG: \" + cg + \" # transversion\\n\" +\n  \"substitution_GT: \" + gt + \" # transversion\\n\" +\n  \"substitution_TG: \" + tg + \" # transversion\\n\")\n}",
            "class": "Expression",
            "engine": "#cwl-js-engine"
          },
          "filename": "config.yaml"
        },
        {
          "fileContent": {
            "script": "{\n  function move_files(files, path){\n    var result = \"\"\n    for(i = 0; i < files.length; i++){\n      result += \"\\nmv \"\n      result += files[i].path + path\n      result += files[i].path.split('/').pop()\n    }\n    result += \"\\n\"\n    return(result)\n  }\n  var result = \"#!/bin/bash\\n\\nmkdir data\"\n  result += move_files([].concat($job.inputs.fastq_files), \" data/\")\n  if($job.inputs.collapsed != null){\n    result += \"mkdir collapsed\"\n    result += move_files([].concat($job.inputs.collapsed), \" collapsed/\")\n  }\n  result += \"source activate quagmir\\nsnakemake -s /opt/quagmir/Snakefile -j \"\n  var nthreads = $job.inputs.number_of_threads || [].concat($job.inputs.fastq_files).length\n  if(nthreads > 36) nthreads = 36\n  result += nthreads\n  \n  return(result)\n}",
            "class": "Expression",
            "engine": "#cwl-js-engine"
          },
          "filename": "run.sh"
        }
      ],
      "class": "CreateFileRequirement"
    },
    {
      "class": "ExpressionEngineRequirement",
      "requirements": [
        {
          "dockerPull": "rabix/js-engine",
          "class": "DockerRequirement"
        }
      ],
      "id": "#cwl-js-engine"
    }
  ],
  "inputs": [
    {
      "type": [
        "null",
        "float"
      ],
      "label": "Minimum ratio",
      "description": "Minimum percentege a sequence comprises out of total number of reads in order for it not to be discarded. Set to -1 in order to disable.",
      "sbg:toolDefaultValue": "-1",
      "id": "#min_ratio",
      "sbg:category": "Filtering parameters"
    },
    {
      "type": [
        "null",
        "int"
      ],
      "label": "Min reads",
      "description": "Minimum coverage a sequence can have in order for it not to be discarded. Set to -1 in order to disable.",
      "sbg:toolDefaultValue": "-1",
      "id": "#min_read",
      "sbg:category": "Filtering parameters"
    },
    {
      "type": [
        "File"
      ],
      "label": "Motif-consensus file",
      "id": "#motif_consensus_file",
      "sbg:fileTypes": "FA",
      "description": "Reference file containing name, motif sequence, and mature miRNA sequence for each miRNA.",
      "sbg:category": "Input files"
    },
    {
      "type": [
        "null",
        "boolean"
      ],
      "label": "Destructive motif pull",
      "description": "When turned on, only the best fitting alignment position is kept, while all other possibilities are discarded.",
      "sbg:toolDefaultValue": "False",
      "id": "#destructive_motif_pull",
      "sbg:category": "Input parameters"
    },
    {
      "type": [
        {
          "items": "File",
          "type": "array"
        }
      ],
      "label": "Input reads",
      "id": "#fastq_files",
      "sbg:fileTypes": "FASTQ",
      "sbg:stageInput": "link",
      "description": "Input miRNA-Seq reads. It can be one sample, or multiple, in which case QuagmiR will use all of available threads so the task's run time would be as short as possible.",
      "sbg:category": "Input files"
    },
    {
      "type": [
        "null",
        "int"
      ],
      "label": "Number of threads",
      "id": "#number_of_threads",
      "description": "Number of threads for QuagmiR to use. If there is 36 files or less, it is recommended to use one thread per file. If there are more files it is best to use divider of the number of files, or in case number of files exceeds 100, set it to 36. Default number of threads is set to number of files in case there is no more than 36 files, and to 36 otherwise.",
      "sbg:category": "Input parameters"
    },
    {
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "label": "Collapsed files",
      "id": "#collapsed",
      "sbg:fileTypes": "COLLAPSED",
      "sbg:stageInput": "link",
      "description": "Collapsed files made by previous tun of the tool. They must be named same like the fastqs with \".collapsed\" at the end of the file name.",
      "sbg:category": "Input files"
    },
    {
      "type": [
        "null",
        "boolean"
      ],
      "label": "Display group output",
      "description": "Option to output group result files for all files.",
      "sbg:toolDefaultValue": "True",
      "id": "#display_group_output",
      "sbg:category": "Output options"
    },
    {
      "type": [
        "null",
        "string"
      ],
      "label": "Group output name",
      "id": "#group_output_name",
      "description": "Prefix for group output files.",
      "sbg:category": "Output options"
    },
    {
      "type": [
        "null",
        "int"
      ],
      "label": "Edit distance filtering treshold on 3' end",
      "description": "Filter out every read that has edit distance on 3' between miRNA sequence and reference greater than set number. To be set to -1 in case filtering on 3' end needs to be swiched off.",
      "sbg:toolDefaultValue": "-1",
      "id": "#edit_distance_treshold_3p",
      "sbg:category": "Filtering parameters"
    },
    {
      "type": [
        "null",
        "int"
      ],
      "label": "Edit distance filtering threshold on 5' end",
      "description": "Filter out every read that has edit distance on 5' between miRNA sequence and reference greater than set number. To be set to -1 in case filtering on 5' end needs to be swiched off.",
      "sbg:toolDefaultValue": "-1",
      "id": "#edit_distance_treshold_5p",
      "sbg:category": "Filtering parameters"
    },
    {
      "type": [
        "null",
        "boolean"
      ],
      "label": "Ambiguous letters support",
      "description": "In case ambiguous letters support is on, the program will be able to work with letters such as N, R, Y, etc, but it will work slower. If it is switched off, all ambiguous letters (even Ns in fastq files) will be treted as missmatches against all bases.",
      "sbg:toolDefaultValue": "False",
      "id": "#ambiguous_letters",
      "sbg:category": "Input parameters"
    },
    {
      "type": [
        "null",
        "float"
      ],
      "label": "Deletion score",
      "description": "Deletion score for calculating edit distance.",
      "sbg:toolDefaultValue": "1",
      "sbg:stageInput": null,
      "id": "#deletion_score",
      "sbg:category": "Edit distance parameters"
    },
    {
      "type": [
        "null",
        "float"
      ],
      "label": "Insertion score",
      "description": "Insertion score for edit distance calculation.",
      "sbg:toolDefaultValue": "1",
      "sbg:stageInput": null,
      "id": "#insertion_score",
      "sbg:category": "Edit distance parameters"
    },
    {
      "type": [
        "null",
        {
          "fields": [
            {
              "type": [
                "null",
                "float"
              ],
              "label": "A->C",
              "name": "ac",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for A to C substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "A->G",
              "name": "ag",
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for A to G substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "A->T",
              "name": "at",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for A to T substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "C->A",
              "name": "ca",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for C to A substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "C->G",
              "name": "cg",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for C to G substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "C->T",
              "name": "ct",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for C to T substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "G->A",
              "name": "ga",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for G to A substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "G->C",
              "name": "gc",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for G to C substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "G->T",
              "name": "gt",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for G to T substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "T->A",
              "name": "ta",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for T to A substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "T->C",
              "name": "tc",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for T to C substitution.",
              "sbg:category": "Edit distance parameters"
            },
            {
              "type": [
                "null",
                "float"
              ],
              "label": "T->G",
              "name": "tg",
              "sbg:toolDefaultValue": "1",
              "sbg:stageInput": null,
              "description": "Edit distance score for T to G substitution.",
              "sbg:category": "Edit distance parameters"
            }
          ],
          "type": "record",
          "name": "substitution_score"
        }
      ],
      "label": "Substitution score",
      "id": "#substitution_score",
      "description": "Substitution score for edit distance calculation.",
      "sbg:category": "Edit distance parameters"
    },
    {
      "type": [
        "File"
      ],
      "label": "Reference file",
      "sbg:fileTypes": "TSV",
      "sbg:stageInput": "link",
      "id": "#reference_file",
      "description": "Reference file."
    },
    {
      "type": [
        "null",
        "string"
      ],
      "label": "miRBase version",
      "id": "#mirbase_version",
      "description": "miRBase version.",
      "sbg:toolDefaultValue": "21"
    },
    {
      "type": [
        "null",
        "string"
      ],
      "label": "Source ontology",
      "id": "#source_ontology",
      "description": "Source ontology.",
      "sbg:toolDefaultValue": "miRBase v21 doi:10.25504/fairsharing.hmgte8"
    }
  ],
  "outputs": [
    {
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "label": "Collapsed fastq files",
      "outputBinding": {
        "glob": {
          "script": "{\n  function file_name(file){\n    return(file.path.split('/').pop())\n  }\n  function get_file_names(input){\n    var input_list = [].concat(input)\n    var input_names = []\n    for(i=0; i<input_list.length; i++){\n      input_names = input_names.concat(file_name(input_list[i]))\n    }\n    return(input_names)\n  }\n  function include(fastq_name, collapsed_names){\n    len = collapsed_names.length\n    var a = 0\n    for(i=0; i< len; i++)\n      if(collapsed_names[i] == fastq_name + \".collapsed\")\n        return(false)\n    return(true)\n  }\n  if($job.inputs.collapsed == null){\n    return(\"collapsed/*\")\n  }\n  else{\n    var result = \"{\"\n    var fastq_names = get_file_names($job.inputs.fastq_files)\n    var collapsed_names = get_file_names($job.inputs.collapsed)\n    for(j=0; j< fastq_names.length; j++){\n      if(include(fastq_names[j], collapsed_names)){\n        result += \"collapsed/\" + fastq_names[j] + \".collapsed,\"\n      }\n    }\n    if(result == \"{\"){\n      return(\"\")\n    }\n    result = result.slice(0, -1) + \"}\"\n    if(result.indexOf(\",\") == -1){\n      result = result.slice(1, -1)\n    }\n    return(result)\n  }\n}",
          "class": "Expression",
          "engine": "#cwl-js-engine"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      },
      "sbg:fileTypes": "COLLAPSED",
      "id": "#collapsed_fastqs",
      "description": "Collapsed fastq files."
    },
    {
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "label": "Sequence info",
      "outputBinding": {
        "glob": "results/*sequence_info.tsv",
        "sbg:inheritMetadataFrom": "#fastq_files"
      },
      "sbg:fileTypes": "TSV",
      "id": "#sequence_info",
      "description": "TSV file containing information about aligned sequences."
    },
    {
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "label": "Nucleotide distribution",
      "outputBinding": {
        "glob": {
          "script": "{\n  if(!($job.inputs.display_nucleotide_dist == false))\n    return(\"results/*nucleotide_dist.tsv\")\n}",
          "class": "Expression",
          "engine": "#cwl-js-engine"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      },
      "sbg:fileTypes": "TSV",
      "id": "#nucleotide_dist",
      "description": "TSV file containing information regarding the distribution of nucleotides."
    },
    {
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "label": "Summary",
      "outputBinding": {
        "glob": {
          "script": "{\n  if(!($job.inputs.display_summary == false))\n    return(\"results/*isomir.tsv\")\n}",
          "class": "Expression",
          "engine": "#cwl-js-engine"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      },
      "sbg:fileTypes": "TSV",
      "id": "#summary",
      "description": "File containing summary information for each reference miRNA sequence (number of reads aligned to it, number of different isomiRs, etc)."
    },
    {
      "type": [
        "null",
        "File"
      ],
      "label": "Group summary",
      "outputBinding": {
        "glob": {
          "script": "{\n  if(!($job.inputs.display_summary == false) &&\n    !($job.inputs.display_group_output == false))\n    return(\"group_results/*isomir.tsv\")\n}",
          "class": "Expression",
          "engine": "#cwl-js-engine"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      },
      "sbg:fileTypes": "TSV",
      "id": "#group_summary",
      "description": "Summary for all inputs."
    },
    {
      "type": [
        "null",
        "File"
      ],
      "label": "Group sequence info",
      "outputBinding": {
        "glob": "group_results/*sequence_info.tsv"
      },
      "sbg:fileTypes": "TSV",
      "id": "#group_sequence_info",
      "description": "Sequence info for all inputs."
    },
    {
      "type": [
        "null",
        "File"
      ],
      "label": "Group nucleotide distribution",
      "outputBinding": {
        "glob": {
          "script": "{\n  if(!($job.inputs.display_nucleotide_dist == false) &&\n    !($job.inputs.display_group_output == false))\n    return(\"group_results/*nucleotide_dist.tsv\")\n}",
          "class": "Expression",
          "engine": "#cwl-js-engine"
        }
      },
      "sbg:fileTypes": "TSV",
      "id": "#group_nucleotide_dist",
      "description": "Nucleotide distribution for all inputs."
    },
    {
      "type": [
        "null",
        {
          "type": "array",
          "items": "File"
        }
      ],
      "label": "GFF",
      "outputBinding": {
        "glob": "results/*gff",
        "sbg:inheritMetadataFrom": "#fastq_files"
      },
      "sbg:fileTypes": "GFF",
      "id": "#gff",
      "description": "GFF."
    }
  ],
  "hints": [
    {
      "value": {
        "script": "{\n  var nthreads = $job.inputs.number_of_threads || [].concat($job.inputs.fastq_files).length\n  if(nthreads > 36) nthreads = 36\n  return(nthreads)\n}",
        "class": "Expression",
        "engine": "#cwl-js-engine"
      },
      "class": "sbg:CPURequirement"
    },
    {
      "value": 1000,
      "class": "sbg:MemRequirement"
    },
    {
      "dockerPull": "images.sbgenomics.com/nikola_tesic/quagmir:1.0",
      "dockerImageId": "",
      "class": "DockerRequirement"
    }
  ],
  "baseCommand": [
    "/bin/bash",
    "run.sh"
  ],
  "stdin": "",
  "stdout": "",
  "successCodes": [],
  "temporaryFailCodes": [],
  "arguments": [],
  "sbg:projectName": "QuagmiR - dev",
  "sbg:toolkit": "QuagmiR",
  "sbg:modifiedOn": 1533498198,
  "sbg:validationErrors": [],
  "sbg:publisher": "sbg",
  "sbg:modifiedBy": "nikola_tesic",
  "sbg:image_url": null,
  "sbg:latestRevision": 53,
  "sbg:categories": [
    "miRNA"
  ],
  "sbg:license": "MIT License",
  "sbg:revisionNotes": "changed docker image to guagmir:1.0",
  "sbg:createdOn": 1516888498,
  "sbg:sbgMaintained": false,
  "sbg:job": {
    "allocatedResources": {
      "cpu": 4,
      "mem": 1000
    },
    "inputs": {
      "mirbase_version": "mirbase_version-string-value",
      "insertion_score": 9.969961420283253,
      "min_ratio": 8.293022266993587,
      "collapsed": [
        {
          "secondaryFiles": [],
          "class": "File",
          "size": 0,
          "path": "/path/to/fastq_file-1.ext.collapsed"
        },
        {
          "secondaryFiles": [],
          "class": "File",
          "size": 0,
          "path": "/path/to/fastq_file-2.ext.not_really"
        },
        {
          "secondaryFiles": [],
          "class": "File",
          "size": 0,
          "path": "/path/to/fastq_file-2.ext.collapsed"
        }
      ],
      "substitution_score": {
        "ac": 7.032371894121999,
        "ca": 5.907251526966097,
        "ct": 10.382560524506955,
        "at": 7.07069141380792,
        "gt": 0.05617499809570303,
        "tg": 9.80306639631734,
        "ga": 10.152746571136447,
        "tc": 5.395801837155783,
        "gc": 5.90324132538586,
        "ta": 5.305421054615811,
        "cg": 8.62896289151937,
        "ag": 9.113533436056201
      },
      "motif_consensus_file": {
        "secondaryFiles": [],
        "class": "File",
        "size": 0,
        "path": "/path/to/motif_consensus_file.ext"
      },
      "destructive_motif_pull": false,
      "group_output_name": "",
      "edit_distance_treshold_5p": 10,
      "edit_distance_treshold_3p": 6,
      "source_ontology": "source_ontology-string-value",
      "number_of_threads": null,
      "min_read": 8,
      "ambiguous_letters": true,
      "fastq_files": [
        {
          "secondaryFiles": [],
          "class": "File",
          "metadata": {
            "sample_group": "test1"
          },
          "size": 0,
          "path": "/path/to/fastq_file-1.ext"
        },
        {
          "secondaryFiles": [],
          "class": "File",
          "metadata": {
            "sample_group": "test2"
          },
          "size": 0,
          "path": "/path/to/fastq_file-2.ext"
        },
        {
          "secondaryFiles": [],
          "class": "File",
          "size": 0,
          "path": "/path/to/fastq_file-3.ext"
        },
        {
          "secondaryFiles": [],
          "class": "File",
          "size": 0,
          "path": "/path/to/fastq_file-4.ext"
        }
      ],
      "reference_file": {
        "size": 0,
        "path": "/path/to/reference_file.ext",
        "class": "File",
        "secondaryFiles": []
      },
      "display_group_output": true,
      "deletion_score": 1.9293315992862135
    }
  },
  "cwlVersion": "sbg:draft-2",
  "sbg:cmdPreview": "/bin/bash run.sh",
  "sbg:links": [
    {
      "label": "QuagmiR Homepage",
      "id": "https://github.com/duxan/quagmir"
    },
    {
      "label": "QuagmiR Source Code",
      "id": "https://github.com/duxan/quagmir"
    },
    {
      "label": "QuagmiR Download",
      "id": "https://github.com/duxan/quagmir"
    },
    {
      "label": "QuagmiR Documentation",
      "id": "https://github.com/duxan/quagmir"
    }
  ],
  "sbg:project": "nikola_tesic/quagmir-dev",
  "sbg:id": "nikola_tesic/quagmir-dev/quagmir/53",
  "sbg:appVersion": [
    "sbg:draft-2"
  ],
  "sbg:createdBy": "nikola_tesic",
  "sbg:revisionsInfo": [
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1516888498,
      "sbg:revision": 0,
      "sbg:revisionNotes": null
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1516888701,
      "sbg:revision": 1,
      "sbg:revisionNotes": "initial version"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1516893139,
      "sbg:revision": 2,
      "sbg:revisionNotes": "beginning of writing of description"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1518486793,
      "sbg:revision": 3,
      "sbg:revisionNotes": null
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519140716,
      "sbg:revision": 4,
      "sbg:revisionNotes": "description update"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519146000,
      "sbg:revision": 5,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519150368,
      "sbg:revision": 6,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519150489,
      "sbg:revision": 7,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519151543,
      "sbg:revision": 8,
      "sbg:revisionNotes": null
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519192255,
      "sbg:revision": 9,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519192536,
      "sbg:revision": 10,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519192702,
      "sbg:revision": 11,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519192782,
      "sbg:revision": 12,
      "sbg:revisionNotes": null
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519192861,
      "sbg:revision": 13,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519192921,
      "sbg:revision": 14,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519193162,
      "sbg:revision": 15,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519193309,
      "sbg:revision": 16,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "uros_sipetic",
      "sbg:modifiedOn": 1519211170,
      "sbg:revision": 17,
      "sbg:revisionNotes": "Small grammar fixes in the description."
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519217650,
      "sbg:revision": 18,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519217919,
      "sbg:revision": 19,
      "sbg:revisionNotes": "added links"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519218142,
      "sbg:revision": 20,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519221133,
      "sbg:revision": 21,
      "sbg:revisionNotes": "description chages"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519223810,
      "sbg:revision": 22,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519223947,
      "sbg:revision": 23,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519224084,
      "sbg:revision": 24,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519224141,
      "sbg:revision": 25,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "uros_sipetic",
      "sbg:modifiedOn": 1519224680,
      "sbg:revision": 26,
      "sbg:revisionNotes": "Small grammar update in description."
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519234204,
      "sbg:revision": 27,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519234267,
      "sbg:revision": 28,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519234341,
      "sbg:revision": 29,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519234374,
      "sbg:revision": 30,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519234580,
      "sbg:revision": 31,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1519236345,
      "sbg:revision": 32,
      "sbg:revisionNotes": "description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533330103,
      "sbg:revision": 33,
      "sbg:revisionNotes": "adding gff test1"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533330244,
      "sbg:revision": 34,
      "sbg:revisionNotes": "changing docker image and adding output file"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533332155,
      "sbg:revision": 35,
      "sbg:revisionNotes": "adding gff test2"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533332465,
      "sbg:revision": 36,
      "sbg:revisionNotes": "adding gff test3"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533333599,
      "sbg:revision": 37,
      "sbg:revisionNotes": "adding gff test4"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533334430,
      "sbg:revision": 38,
      "sbg:revisionNotes": "adding gff test5"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533335019,
      "sbg:revision": 39,
      "sbg:revisionNotes": "adding gff test6"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533335430,
      "sbg:revision": 40,
      "sbg:revisionNotes": "adding gff test7"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533338417,
      "sbg:revision": 41,
      "sbg:revisionNotes": "test8"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533339647,
      "sbg:revision": 42,
      "sbg:revisionNotes": "test9"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533346891,
      "sbg:revision": 43,
      "sbg:revisionNotes": "test10"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533347185,
      "sbg:revision": 44,
      "sbg:revisionNotes": "test11"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533352404,
      "sbg:revision": 45,
      "sbg:revisionNotes": "test12"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533352897,
      "sbg:revision": 46,
      "sbg:revisionNotes": "test13"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533353787,
      "sbg:revision": 47,
      "sbg:revisionNotes": "test14"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533353855,
      "sbg:revision": 48,
      "sbg:revisionNotes": "test15"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533354428,
      "sbg:revision": 49,
      "sbg:revisionNotes": "test16"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533421640,
      "sbg:revision": 50,
      "sbg:revisionNotes": "test17"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533431900,
      "sbg:revision": 51,
      "sbg:revisionNotes": "minor description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533432041,
      "sbg:revision": 52,
      "sbg:revisionNotes": "minor description changes"
    },
    {
      "sbg:modifiedBy": "nikola_tesic",
      "sbg:modifiedOn": 1533498198,
      "sbg:revision": 53,
      "sbg:revisionNotes": "changed docker image to guagmir:1.0"
    }
  ],
  "sbg:revision": 53,
  "sbg:contributors": [
    "nikola_tesic",
    "uros_sipetic"
  ],
  "sbg:toolkitVersion": "1.0",
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "sbg:toolAuthor": "NCI"
}
