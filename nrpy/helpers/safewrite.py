"""
Placeholder

Author: Steven Brandt
"""
from typing import Any, Optional, Union, Dict
import io
import os
from difflib import context_diff
import sys
from subprocess import Popen, PIPE
from nrpy.helpers.colored import colored
from pathlib import Path

try:
    from clang_format import _get_executable as get_executable #type: ignore[import-not-found]
except:
    def get_executable(_:str)->str:
        return "/bin/cat"

clang_formatter = get_executable("clang-format")

verbose = False
nochange = False

class SafeWrite:
    def __init__(self, fname:Union[Path,str], encoding:Optional[str]=None, do_format:bool=False)->None:
        self.fname = os.path.abspath(str(fname))
        self.fd : io.StringIO
        self.do_format = do_format
        self.encoding = encoding
    def __enter__(self)->io.TextIOWrapper:
        self.fd = io.StringIO()
        return self.fd
    def __exit__(self, ty:Any, val:Any, tb:Any)->None:
        print("Checking",self.fname,end="...")
        newcontent = self.fd.getvalue()
        if self.do_format:
            pipe = Popen([clang_formatter],stdout=PIPE,stdin=PIPE,universal_newlines=True)
            out, err = pipe.communicate(newcontent)
            newcontent = out
        if os.path.exists(self.fname):
            kwargs:Dict[str,str] = dict()
            if self.encoding is not None:
                with open(self.fname, encoding=self.encoding) as fd:
                    oldcontent = fd.read()
            else:
                with open(self.fname) as fd:
                    oldcontent = fd.read()
            do_write = newcontent.strip() != oldcontent.strip()
            if do_write and verbose:
                print("Diff for:",self.fname)
                oldlines=[line+"\n" for line in oldcontent.strip().split("\n")]
                newlines=[line+"\n" for line in newcontent.strip().split("\n")]
                sys.stdout.writelines(context_diff(oldlines,newlines,fromfile='before',tofile='after'))
        else:
            do_write = True
        if do_write:
            assert nochange == False
            if self.encoding is not None:
                with open(self.fname, "w", encoding=self.encoding) as fd:
                    fd.write(newcontent)
            else:
                with open(self.fname, "w") as fd:
                    fd.write(newcontent)
            print(colored("[written]","red"))
        else:
            print(colored("[no changes]","green"))
