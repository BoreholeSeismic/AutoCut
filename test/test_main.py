import sys
from autocut.auto_cut import AutoCut
# only for script testing

if __name__ == '__main__':
    autocut = AutoCut()
    autocut.set_input_dir('input')
    autocut.set_output_dir('output')
    autocut.cut()
