import sys
from autocut.auto_cut import AutoCut


if __name__ == '__main__':
    args = sys.argv[1:]
    autocut = AutoCut()
    autocut.set_input_dir('input')
    autocut.set_output_dir('output')
    autocut.cut()
