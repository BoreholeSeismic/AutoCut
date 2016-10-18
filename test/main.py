from auto_cut import AutoCut


if __name__ == '__main__':
    autocut = AutoCut()
    autocut.set_input_dir('input')
    autocut.set_output_dir('output')
    autocut.cut('Cont_Ludwig#9_20160904_021052_553_trace3C.mat')
