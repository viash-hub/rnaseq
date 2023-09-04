workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch

    emit: 
        ouput_ch
}