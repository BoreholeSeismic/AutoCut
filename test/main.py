from auto_cut import AutoCut


if __name__ == '__main__':
    autocut = AutoCut()
    autocut.set_input_dir('input')
    autocut.set_output_dir('output')
    autocut.cut('Cont_Ludwig_20160818_020540_617.mat')
