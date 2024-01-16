from typing_extensions import Literal
from typing import Any, Dict
import sys

color_names = Literal["red","green","yellow","blue","magenta","cyan"]

def not_colored(arg:Any,c:color_names)->str:
    return repr(arg)

colors : Dict[color_names,str] = {
  "red":"\033[31m",
  "green":"\033[32m",
  "yellow":"\033[33m",
  "blue":"\033[34m",
  "magenta":"\033[35m",
  "cyan":"\033[36m",
}
reset = "\033[0m"

def colored(arg:Any,c:color_names)->str:
    assert type(c) == str
    assert c in colors
    s = str(arg)
    return colors[c] + s + reset

if hasattr(sys.stdout,"isatty"):
    is_tty = sys.stdout.isatty()
else:
    is_tty = False

is_jupyter = type(sys.stdout).__name__ == 'OutStream' and  type(sys.stdout).__module__ == 'ipykernel.iostream'
if (not is_tty) and (not is_jupyter):
    colored = not_colored

if __name__ == "__main__":
    print(colored("red", "red"), colored("green", "green"))
