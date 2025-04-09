rule move_files_to_results:
    input:
        summary = "resources/table_{version}/" + RUN_NAME + ".countsummary.txt",
    output:
        summary = "resources/table_{version}/results/" + RUN_NAME + ".countsummary.txt",
    shadow: "shallow"
    shell: 
        "mv {input.summary} {output.summary}"
