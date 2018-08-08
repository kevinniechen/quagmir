{
  "stdin": "",
  "outputs": [
    {
      "label": "Collapsed fastq files",
      "sbg:fileTypes": "COLLAPSED",
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "description": "Collapsed fastq files.",
      "id": "#collapsed_fastqs",
      "outputBinding": {
        "glob": {
          "engine": "#cwl-js-engine",
          "class": "Expression",
          "script": "{\n  function file_name(file){\n    return(file.path.split('/').pop())\n  }\n  function get_file_names(input){\n    var input_list = [].concat(input)\n    var input_names = []\n    for(i=0; i<input_list.length; i++){\n      input_names = input_names.concat(file_name(input_list[i]))\n    }\n    return(input_names)\n  }\n  function include(fastq_name, collapsed_names){\n    len = collapsed_names.length\n    var a = 0\n    for(i=0; i< len; i++)\n      if(collapsed_names[i] == fastq_name + \".collapsed\")\n        return(false)\n    return(true)\n  }\n  if($job.inputs.collapsed == null){\n    return(\"collapsed/*\")\n  }\n  else{\n    var result = \"{\"\n    var fastq_names = get_file_names($job.inputs.fastq_files)\n    var collapsed_names = get_file_names($job.inputs.collapsed)\n    for(j=0; j< fastq_names.length; j++){\n      if(include(fastq_names[j], collapsed_names)){\n        result += \"collapsed/\" + fastq_names[j] + \".collapsed,\"\n      }\n    }\n    if(result == \"{\"){\n      return(\"\")\n    }\n    result = result.slice(0, -1) + \"}\"\n    if(result.indexOf(\",\") == -1){\n      result = result.slice(1, -1)\n    }\n    return(result)\n  }\n}"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      }
    },
    {
      "label": "Sequence info",
      "sbg:fileTypes": "TSV",
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "description": "TSV file containing information about aligned sequences.",
      "id": "#sequence_info",
      "outputBinding": {
        "glob": "results/*sequence_info.tsv",
        "sbg:inheritMetadataFrom": "#fastq_files"
      }
    },
    {
      "label": "Nucleotide distribution",
      "sbg:fileTypes": "TSV",
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "description": "TSV file containing information regarding the distribution of nucleotides.",
      "id": "#nucleotide_dist",
      "outputBinding": {
        "glob": {
          "engine": "#cwl-js-engine",
          "class": "Expression",
          "script": "{\n  if(!($job.inputs.display_nucleotide_dist == false))\n    return(\"results/*nucleotide_dist.tsv\")\n}"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      }
    },
    {
      "label": "Summary",
      "sbg:fileTypes": "TSV",
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "description": "File containing summary information for each reference miRNA sequence (number of reads aligned to it, number of different isomiRs, etc).",
      "id": "#summary",
      "outputBinding": {
        "glob": {
          "engine": "#cwl-js-engine",
          "class": "Expression",
          "script": "{\n  if(!($job.inputs.display_summary == false))\n    return(\"results/*isomir.tsv\")\n}"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      }
    },
    {
      "label": "Group summary",
      "sbg:fileTypes": "TSV",
      "type": [
        "null",
        "File"
      ],
      "description": "Summary for all inputs.",
      "id": "#group_summary",
      "outputBinding": {
        "glob": {
          "engine": "#cwl-js-engine",
          "class": "Expression",
          "script": "{\n  if(!($job.inputs.display_summary == false) &&\n    !($job.inputs.display_group_output == false))\n    return(\"group_results/*isomir.tsv\")\n}"
        },
        "sbg:inheritMetadataFrom": "#fastq_files"
      }
    },
    {
      "label": "Group sequence info",
      "sbg:fileTypes": "TSV",
      "type": [
        "null",
        "File"
      ],
      "description": "Sequence info for all inputs.",
      "id": "#group_sequence_info",
      "outputBinding": {
        "glob": "group_results/*sequence_info.tsv"
      }
    },
    {
      "label": "Group nucleotide distribution",
      "sbg:fileTypes": "TSV",
      "type": [
        "null",
        "File"
      ],
      "description": "Nucleotide distribution for all inputs.",
      "id": "#group_nucleotide_dist",
      "outputBinding": {
        "glob": {
          "engine": "#cwl-js-engine",
          "class": "Expression",
          "script": "{\n  if(!($job.inputs.display_nucleotide_dist == false) &&\n    !($job.inputs.display_group_output == false))\n    return(\"group_results/*nucleotide_dist.tsv\")\n}"
        }
      }
    },
    {
      "label": "GFF",
      "sbg:fileTypes": "GFF",
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "description": "GFF.",
      "id": "#gff",
      "outputBinding": {
        "glob": "results/*gff",
        "sbg:inheritMetadataFrom": "#fastq_files"
      }
    }
  ],
  "sbg:toolkitVersion": "1.0",
  "inputs": [
    {
      "label": "Minimum ratio",
      "sbg:toolDefaultValue": "-1",
      "id": "#min_ratio",
      "description": "Minimum percentege a sequence comprises out of total number of reads in order for it not to be discarded. Set to -1 in order to disable.",
      "type": [
        "null",
        "float"
      ],
      "sbg:category": "Filtering parameters"
    },
    {
      "label": "Min reads",
      "sbg:toolDefaultValue": "-1",
      "id": "#min_read",
      "description": "Minimum coverage a sequence can have in order for it not to be discarded. Set to -1 in order to disable.",
      "type": [
        "null",
        "int"
      ],
      "sbg:category": "Filtering parameters"
    },
    {
      "label": "Motif-consensus file",
      "sbg:fileTypes": "FA",
      "id": "#motif_consensus_file",
      "description": "Reference file containing name, motif sequence, and mature miRNA sequence for each miRNA.",
      "type": [
        "File"
      ],
      "sbg:category": "Input files"
    },
    {
      "label": "Destructive motif pull",
      "sbg:toolDefaultValue": "False",
      "id": "#destructive_motif_pull",
      "description": "When turned on, only the best fitting alignment position is kept, while all other possibilities are discarded.",
      "type": [
        "null",
        "boolean"
      ],
      "sbg:category": "Input parameters"
    },
    {
      "label": "Input reads",
      "sbg:stageInput": "link",
      "sbg:fileTypes": "FASTQ",
      "id": "#fastq_files",
      "description": "Input miRNA-Seq reads. It can be one sample, or multiple, in which case QuagmiR will use all of available threads so the task's run time would be as short as possible.",
      "type": [
        {
          "items": "File",
          "type": "array"
        }
      ],
      "sbg:category": "Input files"
    },
    {
      "label": "Number of threads",
      "description": "Number of threads for QuagmiR to use. If there is 36 files or less, it is recommended to use one thread per file. If there are more files it is best to use divider of the number of files, or in case number of files exceeds 100, set it to 36. Default number of threads is set to number of files in case there is no more than 36 files, and to 36 otherwise.",
      "type": [
        "null",
        "int"
      ],
      "id": "#number_of_threads",
      "sbg:category": "Input parameters"
    },
    {
      "label": "Collapsed files",
      "sbg:stageInput": "link",
      "sbg:fileTypes": "COLLAPSED",
      "id": "#collapsed",
      "description": "Collapsed files made by previous tun of the tool. They must be named same like the fastqs with \".collapsed\" at the end of the file name.",
      "type": [
        "null",
        {
          "items": "File",
          "type": "array"
        }
      ],
      "sbg:category": "Input files"
    },
    {
      "label": "Display group output",
      "sbg:toolDefaultValue": "True",
      "id": "#display_group_output",
      "description": "Option to output group result files for all files.",
      "type": [
        "null",
        "boolean"
      ],
      "sbg:category": "Output options"
    },
    {
      "label": "Group output name",
      "description": "Prefix for group output files.",
      "type": [
        "null",
        "string"
      ],
      "id": "#group_output_name",
      "sbg:category": "Output options"
    },
    {
      "label": "Edit distance filtering treshold on 3' end",
      "sbg:toolDefaultValue": "-1",
      "id": "#edit_distance_treshold_3p",
      "description": "Filter out every read that has edit distance on 3' between miRNA sequence and reference greater than set number. To be set to -1 in case filtering on 3' end needs to be swiched off.",
      "type": [
        "null",
        "int"
      ],
      "sbg:category": "Filtering parameters"
    },
    {
      "label": "Edit distance filtering threshold on 5' end",
      "sbg:toolDefaultValue": "-1",
      "id": "#edit_distance_treshold_5p",
      "description": "Filter out every read that has edit distance on 5' between miRNA sequence and reference greater than set number. To be set to -1 in case filtering on 5' end needs to be swiched off.",
      "type": [
        "null",
        "int"
      ],
      "sbg:category": "Filtering parameters"
    },
    {
      "label": "Ambiguous letters support",
      "sbg:toolDefaultValue": "False",
      "id": "#ambiguous_letters",
      "description": "In case ambiguous letters support is on, the program will be able to work with letters such as N, R, Y, etc, but it will work slower. If it is switched off, all ambiguous letters (even Ns in fastq files) will be treted as missmatches against all bases.",
      "type": [
        "null",
        "boolean"
      ],
      "sbg:category": "Input parameters"
    },
    {
      "label": "Deletion score",
      "sbg:stageInput": null,
      "sbg:toolDefaultValue": "1",
      "id": "#deletion_score",
      "description": "Deletion score for calculating edit distance.",
      "type": [
        "null",
        "float"
      ],
      "sbg:category": "Edit distance parameters"
    },
    {
      "label": "Insertion score",
      "sbg:stageInput": null,
      "sbg:toolDefaultValue": "1",
      "id": "#insertion_score",
      "description": "Insertion score for edit distance calculation.",
      "type": [
        "null",
        "float"
      ],
      "sbg:category": "Edit distance parameters"
    },
    {
      "label": "Substitution score",
      "description": "Substitution score for edit distance calculation.",
      "type": [
        "null",
        {
          "name": "substitution_score",
          "type": "record",
          "fields": [
            {
              "label": "A->C",
              "name": "ac",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for A to C substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "A->G",
              "name": "ag",
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for A to G substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "A->T",
              "name": "at",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for A to T substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "C->A",
              "name": "ca",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for C to A substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "C->G",
              "name": "cg",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for C to G substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "C->T",
              "name": "ct",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for C to T substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "G->A",
              "name": "ga",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for G to A substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "G->C",
              "name": "gc",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for G to C substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "G->T",
              "name": "gt",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for G to T substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "T->A",
              "name": "ta",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for T to A substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "T->C",
              "name": "tc",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for T to C substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            },
            {
              "label": "T->G",
              "name": "tg",
              "sbg:stageInput": null,
              "sbg:toolDefaultValue": "1",
              "description": "Edit distance score for T to G substitution.",
              "type": [
                "null",
                "float"
              ],
              "sbg:category": "Edit distance parameters"
            }
          ]
        }
      ],
      "id": "#substitution_score",
      "sbg:category": "Edit distance parameters"
    },
    {
      "label": "Reference file",
      "sbg:stageInput": "link",
      "sbg:fileTypes": "TSV",
      "id": "#reference_file",
      "description": "Reference file.",
      "type": [
        "File"
      ]
    },
    {
      "label": "miRBase version",
      "description": "miRBase version.",
      "type": [
        "null",
        "string"
      ],
      "id": "#mirbase_version",
      "sbg:toolDefaultValue": "21"
    },
    {
      "label": "Source ontology",
      "description": "Source ontology.",
      "type": [
        "null",
        "string"
      ],
      "id": "#source_ontology",
      "sbg:toolDefaultValue": "miRBase v21 doi:10.25504/fairsharing.hmgte8"
    }
  ],
  "sbg:job": {
    "allocatedResources": {
      "cpu": 4,
      "mem": 1000
    },
    "inputs": {
      "source_ontology": "source_ontology-string-value",
      "reference_file": {
        "path": "/path/to/reference_file.ext",
        "size": 0,
        "secondaryFiles": [],
        "class": "File"
      },
      "number_of_threads": null,
      "display_group_output": true,
      "insertion_score": 9.969961420283253,
      "fastq_files": [
        {
          "path": "/path/to/fastq_file-1.ext",
          "size": 0,
          "metadata": {
            "sample_group": "test1"
          },
          "secondaryFiles": [],
          "class": "File"
        },
        {
          "path": "/path/to/fastq_file-2.ext",
          "size": 0,
          "metadata": {
            "sample_group": "test2"
          },
          "secondaryFiles": [],
          "class": "File"
        },
        {
          "path": "/path/to/fastq_file-3.ext",
          "size": 0,
          "secondaryFiles": [],
          "class": "File"
        },
        {
          "path": "/path/to/fastq_file-4.ext",
          "size": 0,
          "secondaryFiles": [],
          "class": "File"
        }
      ],
      "min_read": 8,
      "edit_distance_treshold_3p": 6,
      "mirbase_version": "mirbase_version-string-value",
      "min_ratio": 8.293022266993587,
      "substitution_score": {
        "tc": 5.395801837155783,
        "ct": 10.382560524506955,
        "cg": 8.62896289151937,
        "at": 7.07069141380792,
        "ac": 7.032371894121999,
        "ga": 10.152746571136447,
        "gt": 0.05617499809570303,
        "gc": 5.90324132538586,
        "ag": 9.113533436056201,
        "tg": 9.80306639631734,
        "ta": 5.305421054615811,
        "ca": 5.907251526966097
      },
      "collapsed": [
        {
          "path": "/path/to/fastq_file-1.ext.collapsed",
          "size": 0,
          "secondaryFiles": [],
          "class": "File"
        },
        {
          "path": "/path/to/fastq_file-2.ext.not_really",
          "size": 0,
          "secondaryFiles": [],
          "class": "File"
        },
        {
          "path": "/path/to/fastq_file-2.ext.collapsed",
          "size": 0,
          "secondaryFiles": [],
          "class": "File"
        }
      ],
      "edit_distance_treshold_5p": 10,
      "motif_consensus_file": {
        "path": "/path/to/motif_consensus_file.ext",
        "size": 0,
        "secondaryFiles": [],
        "class": "File"
      },
      "group_output_name": "",
      "deletion_score": 1.9293315992862135,
      "ambiguous_letters": true,
      "destructive_motif_pull": false
    }
  },
  "sbg:toolkit": "QuagmiR",
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
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "hints": [
    {
      "value": {
        "engine": "#cwl-js-engine",
        "class": "Expression",
        "script": "{\n  var nthreads = $job.inputs.number_of_threads || [].concat($job.inputs.fastq_files).length\n  if(nthreads > 36) nthreads = 36\n  return(nthreads)\n}"
      },
      "class": "sbg:CPURequirement"
    },
    {
      "value": 1000,
      "class": "sbg:MemRequirement"
    },
    {
      "dockerImageId": "",
      "dockerPull": "images.sbgenomics.com/nikola_tesic/quagmir:1.0",
      "class": "DockerRequirement"
    }
  ],
  "baseCommand": [
    "/bin/bash",
    "run.sh"
  ],
  "id": "https://api.sbgenomics.com/v2/apps/nikola_tesic/quagmir-dev/quagmir/54/raw/",
  "sbg:modifiedBy": "nikola_tesic",
  "temporaryFailCodes": [],
  "sbg:toolAuthor": "NCI",
  "requirements": [
    {
      "fileDef": [
        {
          "filename": "config.yaml",
          "fileContent": {
            "engine": "#cwl-js-engine",
            "class": "Expression",
            "script": "{\n  var motif_consensus_file = $job.inputs.motif_consensus_file.path\n  var reference_file = $job.inputs.reference_file.path\n  var min_ratio = $job.inputs.min_ratio || -1\n  var min_read = $job.inputs.min_read || -1\n  var destructive_motif_pull = $job.inputs.destructive_motif_pull || \"False\"\n  if(destructive_motif_pull == true)destructive_motif_pull = \"True\"\n  var ambiguous_letters = $job.inputs.ambiguous_letters || \"False\"\n  if(ambiguous_letters == true)ambiguous_letters = \"True\"\n  var display_summary = $job.inputs.display_summary || \"False\"\n  if($job.inputs.display_summary == false) var display_summary = \"False\"\n  else var display_summary = \"True\"\n  if($job.inputs.display_sequence_info == false) var display_sequence_info = \"False\"\n  else var display_sequence_info = \"True\"  \n  if($job.inputs.display_nucleotide_dist == false) var display_nucleotide_dist = \"False\"\n  else var display_nucleotide_dist = \"True\"\n  if($job.inputs.display_group_output == false) var display_group_output = \"False\"\n  else var display_group_output = \"True\"\n  if($job.inputs.group_output_name == null ||\n     $job.inputs.group_output_name == \"\") var group_output_name = \"group_output\"\n  else var group_output_name = $job.inputs.group_output_name\n  var mirbase_version = $job.inputs.mirbase_version || \"21\"\n  var source_ontology = $job.inputs.source_ontology || \"miRBase v21 doi:10.25504/fairsharing.hmgte8\"\n  \n  var edit_distance_3p = $job.inputs.edit_distance_treshold_3p || -1\n  var edit_distance_5p = $job.inputs.edit_distance_treshold_5p || -1\n  \n  //edit distance parameters\n  var deletion_score = $job.inputs.deletion_score || 1\n  var insertion_score = $job.inputs.insertion_score || 1\n  if($job.inputs.substitution_score){\n    var ac = $job.inputs.substitution_score.ac || 1\n    var ag = $job.inputs.substitution_score.ag || 1\n    var at = $job.inputs.substitution_score.at || 1\n    var ca = $job.inputs.substitution_score.ca || 1\n    var cg = $job.inputs.substitution_score.cg || 1\n    var ct = $job.inputs.substitution_score.ct || 1\n    var ga = $job.inputs.substitution_score.ga || 1\n    var gc = $job.inputs.substitution_score.gc || 1\n    var gt = $job.inputs.substitution_score.gt || 1\n    var ta = $job.inputs.substitution_score.ta || 1\n    var tc = $job.inputs.substitution_score.tc || 1\n    var tg = $job.inputs.substitution_score.tg || 1\n  }\n  else{\n    var ac = 1\n    var ag = 1\n    var at = 1\n    var ca = 1\n    var cg = 1\n    var ct = 1\n    var ga = 1\n    var gc = 1\n    var gt = 1\n    var ta = 1\n    var tc = 1\n    var tg = 1\n  }\n  \n  return(\"# FILTERING\\n\" + \n  \"min_ratio: \" + min_ratio + \"\\n\" + //.1 \n  \"min_read: \" + min_read + \"\\n\" + //9\n  \"edit_distance_3p: \" + edit_distance_3p + \"\\n\" + //-1\n  \"edit_distance_5p: \" + edit_distance_5p + \"\\n\\n\" + //-1\n\n  \"# DISPLAY\\n\" +\n  \"display_summary: \" + display_summary + \"\\n\" +//True\n  \"display_sequence_info: \" + display_sequence_info + \"\\n\" + //True\n  \"display_nucleotide_dist: \" + display_nucleotide_dist + \"\\n\" + //True\n  \"display_group_output: \" + display_group_output + \"\\n\\n\" + //True\n         \n  \"# AMBIGUOUS LETTERS\\n\" +\n  \"ambiguous_letters: \" + ambiguous_letters + \"\\n\\n\" + //True\n         \n  \"# DESTRUCTIVE MOTIF PULL\\n\" +\n  \"destructive_motif_pull: \" + destructive_motif_pull + \"\\n\\n\" + //False\n\n  \"# INPUT\\n\" +\n  \"data: data\\n\" + // input folder\n  \"motif_consensus_file: \" + motif_consensus_file + \"\\n\" +\n  \"reference_file: \" + reference_file + '\\n' +\n  \"mirbase_version: \" + mirbase_version + \"\\n\" + // version of miRBase being used for reference files\n  \"source_ontology: \" + source_ontology + \"\\n\\n\" +\n  \n  \"# OUTPUT\\n\" +\n  \"group_output_name: \" + group_output_name + \"\\n\\n\" +\n  \n  \"# DISTANCE METRIC\\n\" +\n  \"deletion_score: \" + deletion_score + \"\\n\" +\n  \"insertion_score: \" + insertion_score + \"\\n\" +\n  \"substitution_AG: \" + ag + \" # transition\\n\" +\n  \"substitution_GA: \" + ga + \" # transition\\n\" +\n  \"substitution_CT: \" + ct + \" # transition\\n\" +\n  \"substitution_TC: \" + tc + \" # transition\\n\" +\n  \"substitution_AT: \" + at + \" # transversion\\n\" +\n  \"substitution_TA: \" + ta + \" # transversion\\n\" +\n  \"substitution_AC: \" + ac + \" # transversion\\n\" +\n  \"substitution_CA: \" + ca + \" # transversion\\n\" +\n  \"substitution_GC: \" + gc + \" # transversion\\n\" +\n  \"substitution_CG: \" + cg + \" # transversion\\n\" +\n  \"substitution_GT: \" + gt + \" # transversion\\n\" +\n  \"substitution_TG: \" + tg + \" # transversion\\n\")\n}"
          }
        },
        {
          "filename": "run.sh",
          "fileContent": {
            "engine": "#cwl-js-engine",
            "class": "Expression",
            "script": "{\n  function move_files(files, path){\n    var result = \"\"\n    for(i = 0; i < files.length; i++){\n      result += \"\\nmv \"\n      result += files[i].path + path\n      result += files[i].path.split('/').pop()\n    }\n    result += \"\\n\"\n    return(result)\n  }\n  var result = \"#!/bin/bash\\n\\nmkdir data\"\n  result += move_files([].concat($job.inputs.fastq_files), \" data/\")\n  if($job.inputs.collapsed != null){\n    result += \"mkdir collapsed\"\n    result += move_files([].concat($job.inputs.collapsed), \" collapsed/\")\n  }\n  result += \"source activate quagmir\\nsnakemake -s /opt/quagmir/Snakefile -j \"\n  var nthreads = $job.inputs.number_of_threads || [].concat($job.inputs.fastq_files).length\n  if(nthreads > 36) nthreads = 36\n  result += nthreads\n  \n  return(result)\n}"
          }
        }
      ],
      "class": "CreateFileRequirement"
    },
    {
      "requirements": [
        {
          "dockerPull": "rabix/js-engine",
          "class": "DockerRequirement"
        }
      ],
      "id": "#cwl-js-engine",
      "class": "ExpressionEngineRequirement"
    }
  ],
  "sbg:categories": [
    "miRNA"
  ],
  "class": "CommandLineTool",
  "sbg:project": "nikola_tesic/quagmir-dev",
  "sbg:validationErrors": [],
  "sbg:sbgMaintained": false,
  "label": "QuagmiR",
  "successCodes": [],
  "sbg:license": "MIT License",
  "sbg:createdOn": 1516888498,
  "sbg:publisher": "sbg",
  "sbg:createdBy": "nikola_tesic",
  "sbg:id": "nikola_tesic/quagmir-dev/quagmir/54",
  "sbg:latestRevision": 54,
  "sbg:appVersion": [
    "sbg:draft-2"
  ],
  "sbg:cmdPreview": "/bin/bash run.sh",
  "sbg:projectName": "QuagmiR - dev",
  "sbg:modifiedOn": 1533730848,
  "sbg:revisionsInfo": [
    {
      "sbg:modifiedOn": 1516888498,
      "sbg:revisionNotes": null,
      "sbg:revision": 0,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1516888701,
      "sbg:revisionNotes": "initial version",
      "sbg:revision": 1,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1516893139,
      "sbg:revisionNotes": "beginning of writing of description",
      "sbg:revision": 2,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1518486793,
      "sbg:revisionNotes": null,
      "sbg:revision": 3,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519140716,
      "sbg:revisionNotes": "description update",
      "sbg:revision": 4,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519146000,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 5,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519150368,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 6,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519150489,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 7,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519151543,
      "sbg:revisionNotes": null,
      "sbg:revision": 8,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519192255,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 9,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519192536,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 10,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519192702,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 11,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519192782,
      "sbg:revisionNotes": null,
      "sbg:revision": 12,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519192861,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 13,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519192921,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 14,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519193162,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 15,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519193309,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 16,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519211170,
      "sbg:revisionNotes": "Small grammar fixes in the description.",
      "sbg:revision": 17,
      "sbg:modifiedBy": "uros_sipetic"
    },
    {
      "sbg:modifiedOn": 1519217650,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 18,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519217919,
      "sbg:revisionNotes": "added links",
      "sbg:revision": 19,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519218142,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 20,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519221133,
      "sbg:revisionNotes": "description chages",
      "sbg:revision": 21,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519223810,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 22,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519223947,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 23,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519224084,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 24,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519224141,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 25,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519224680,
      "sbg:revisionNotes": "Small grammar update in description.",
      "sbg:revision": 26,
      "sbg:modifiedBy": "uros_sipetic"
    },
    {
      "sbg:modifiedOn": 1519234204,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 27,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519234267,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 28,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519234341,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 29,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519234374,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 30,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519234580,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 31,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1519236345,
      "sbg:revisionNotes": "description changes",
      "sbg:revision": 32,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533330103,
      "sbg:revisionNotes": "adding gff test1",
      "sbg:revision": 33,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533330244,
      "sbg:revisionNotes": "changing docker image and adding output file",
      "sbg:revision": 34,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533332155,
      "sbg:revisionNotes": "adding gff test2",
      "sbg:revision": 35,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533332465,
      "sbg:revisionNotes": "adding gff test3",
      "sbg:revision": 36,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533333599,
      "sbg:revisionNotes": "adding gff test4",
      "sbg:revision": 37,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533334430,
      "sbg:revisionNotes": "adding gff test5",
      "sbg:revision": 38,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533335019,
      "sbg:revisionNotes": "adding gff test6",
      "sbg:revision": 39,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533335430,
      "sbg:revisionNotes": "adding gff test7",
      "sbg:revision": 40,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533338417,
      "sbg:revisionNotes": "test8",
      "sbg:revision": 41,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533339647,
      "sbg:revisionNotes": "test9",
      "sbg:revision": 42,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533346891,
      "sbg:revisionNotes": "test10",
      "sbg:revision": 43,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533347185,
      "sbg:revisionNotes": "test11",
      "sbg:revision": 44,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533352404,
      "sbg:revisionNotes": "test12",
      "sbg:revision": 45,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533352897,
      "sbg:revisionNotes": "test13",
      "sbg:revision": 46,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533353787,
      "sbg:revisionNotes": "test14",
      "sbg:revision": 47,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533353855,
      "sbg:revisionNotes": "test15",
      "sbg:revision": 48,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533354428,
      "sbg:revisionNotes": "test16",
      "sbg:revision": 49,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533421640,
      "sbg:revisionNotes": "test17",
      "sbg:revision": 50,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533431900,
      "sbg:revisionNotes": "minor description changes",
      "sbg:revision": 51,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533432041,
      "sbg:revisionNotes": "minor description changes",
      "sbg:revision": 52,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533498198,
      "sbg:revisionNotes": "changed docker image to guagmir:1.0",
      "sbg:revision": 53,
      "sbg:modifiedBy": "nikola_tesic"
    },
    {
      "sbg:modifiedOn": 1533730848,
      "sbg:revisionNotes": "minor description changes",
      "sbg:revision": 54,
      "sbg:modifiedBy": "nikola_tesic"
    }
  ],
  "sbg:revision": 54,
  "sbg:contributors": [
    "uros_sipetic",
    "nikola_tesic"
  ],
  "abg:revisionNotes": "minor description changes",
  "arguments": [],
  "stdout": "",
  "description": "**QuagmiR** is a python tool for mapping of miRNA sequences and detection and quantification of different isomiRs. For mapping, the tool searches for exact match in the conserved middle part of the miRNA sequence, and once the match is found it uses filtering based on other information to decide whether to discard the sequence or not. It uses FASTQ files with miRNA-Seq reads as input.\n\n**QuagmiR** requires the following input files and parameters:\n\n* **Input reads**: MiRNA-Seq reads in FASTQ format (gzipped input accepted). If multiple files are provided, all available threads on the instance will be used so the task's run time is made as short as possible by running **QuagmiR** on samples in parallel.\n\n* **Motif-consensus file**: Reference file containing names of miRNA sequences, their conserved parts and miRNA sequences themselves. For this, [motif_list_hsa.fa](https://cgc.sbgenomics.com/public/files/5a8dab554f0c3a81d19fca69/) can be used, which contains human miRNAs. Otherwise, any list of miRNAs can be used as **Reference file** as long as it adheres to the format of [motif_list_hsa.fa](https://cgc.sbgenomics.com/public/files/5a8dab554f0c3a81d19fca69/) (for example, if there is interest only in specific miRNAs, only those can be used as reference, and if FASTQ files belong to some non-human species, then reference file for that species is needed).\n\n* **Reference file**: Reference TSV containing additional information about miRNA references. Example file can be found on QuagmiR's github repo.\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the end of the page.*\n\n### Common Use Cases\n\nMicroRNAs (abbreviated miRNAs) are small non-coding RNA molecules (containing about 22 nucleotides) found in plants, animals and some viruses, that have a function in RNA silencing and post-transcriptional regulation of gene expression. IsomiRs are sequences that have variations with respect to the reference MiRNA sequence.\n\nEach miRNA is considered to have a conserved part in the middle that is the same for all of its isomiRs. **QuagmiR** searches for that part in the input FASTQ file, and than chooses reads among those that pass the filtering thresholds (such as 3' and 5' ends filtering where the ends need to have smaller weighted Levenstein distance from the miRNA than the user defined threshold (which can be switched off altogether)).\n\nAfter filtering, for sequences that are aligned to multiple miRNA reads, the best alignment is chosen based on the total weighted Levenstein distance between them. This can be changed by setting the **Destructive motif pull** parameter to False. That way, all possible alignments will be outputted.\n\nAdditionally, as previously mentioned, different **Reference file**s can be chosen, so that they better fit a specific use-case (i.e. a subset of miRNA references contained in the **Reference file** we provide on the platform, so that only particular types of miRNA that are of interest are explored).\n\nAll these parameters, and several others that are described at the end of this page, provide the user with a wide range of possibilities to use this tool in. The tool can be used for general miRNA mapping and isomiR quantification, but it can be used for more specific needs as well by applying different filters, using different sets of miRNA references, turning the **Destructive motif pull** option on or off, etc.\n\n### Changes Introduced by Seven Bridges\n\n* Collapsed files can be used as an optional input, so for each FASTQ file that has it's COLLAPSED file, **QuagmiR** will skip the collapsing step, so the execution time can be shortened. This is useful when analyzing same input files with different tool options such as filtering thresholds, so the files are collapsed only the first time. If some of the input files have COLLAPSED files and some don't, COLLAPSED files will be created and outputed for the new input files.\n\n### Common Issues and Important Notes\n\n* At the moment, only **Reference file**s for human miRNA have been built and tested, so the tool can't be used on other species without first creating the **Reference file**, which requires information regarding that species' miRNA that may not be known yet (most notably, length and positions of conserved parts within the known miRNA sequences).\n\n### Performance Benchmarking\n\nIn the following table you can find estimates of running times and costs. All input files are FASTQ files containing miRNA-Seq reads. If not specified otherwise, number of threads for the tool will be set to be the same as number of input files, or 36 in case the number of files is greater than 36. The instance is set to either default c4.2xlarge, or, in case more threads are required, to the cheapest instance that has sufficient number of threads.\n\nExecution time depends on many factors, but the most important ones are size of input files and whether the **Ambiguous letters support** is turned on or not.\n\n*Cost can be significantly reduced by **spot instance** usage. Visit [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*         \n\n#####Example runs:\n\n| Experiment type | Input size | Ambiguous letters | Read length | Duration | Cost | Instance (AWS)| Instance type |\n|-----------------------|------------|-----------------|------------|-----------------|-------------|--------------|------------------|-------------|---------------|------------|\n| miRNA-Seq         | 1    | 864 MB     | No        | 15-30         | 17m   | $0.16            | c4.2xlarge      | On-demand        |\n| miRNA-Seq         | 1    | 864 MB     | Yes        | 15-30         | 1h 19m   | $0.73            | c4.2xlarge      | On-demand        |\n| miRNA-Seq         | 72     | 599MB - 872MB      | Yes        | 15-30         | 2h 19m   | $1.24            | c4.8xlarge      | On-demand        |\n\n### API Python Implementation\n\nThe tool's draft task can also be submitted via the **API**. In order to learn how to get your **Authentication token** and **API endpoint** for corresponding platform visit our [documentation](https://github.com/sbg/sevenbridges-python#authentication-and-configuration).\n\n```python\nfrom sevenbridges import Api\n\nauthentication_token, api_endpoint = \"enter_your_token\", \"enter_api_endpoint\"\napi = Api(token=authentication_token, url=api_endpoint)\n# Get project_id/workflow_id from your address bar. Example: https://igor.sbgenomics.com/u/your_username/project/tool\nproject_id, tool_id = \"your_username/project\", \"your_username/project/tool\"\n# Get file names from files in your project. Example: Names are taken from Data/Public Reference Files.\ninputs = {\n    'fastq_files': list(api.files.query(project=project_id, names=['miRNA_test_file.fastq'])), \n    'motif_consensus_file': api.files.query(project=project_id, names=['motif_list_hsa.fa'])[0],\n    'reference_file': api.files.query(project=project_id, names=['miRBase21-master.tsv'])[0]\n}\ntask = api.tasks.create(name='QuagmiR - API Example', project=project_id, app=tool_id, inputs=inputs, run=True)\n```\n\nInstructions for installing and configuring the API Python client, are provided on [github](https://github.com/sbg/sevenbridges-python#installation). For more information about using the API Python client, consult [sevenbridges-python documentation](http://sevenbridges-python.readthedocs.io/en/latest/). **More examples** are available [here](https://github.com/sbg/okAPI).\n\nAdditionally, [API R](https://github.com/sbg/sevenbridges-r) and [API Java](https://github.com/sbg/sevenbridges-java) clients are available. To learn more about using these API clients please refer to the [API R client documentation](https://sbg.github.io/sevenbridges-r/), and [API Java client documentation](https://docs.sevenbridges.com/docs/java-library-quickstart).\n\n### References\n\n[1] [Github repo](https://github.com/duxan/quagmir)\n\n[2] Authors and their contact: [Xavier Bofill De Ross](https://github.com/xbdr86), [Kevin Chen](https://github.com/kevchn), [Nikola Tesic](https://github.com/nikola-tesic), [Dusan Randjelovic](https://github.com/duxan)",
  "sbg:image_url": null,
  "cwlVersion": "sbg:draft-2"
}
