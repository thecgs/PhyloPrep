import argparse, json, runpy, sys
from pathlib import Path
class Captured(BaseException):
    def __init__(self, parser): self.parser=parser
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
original=argparse.ArgumentParser.parse_args
argparse.ArgumentParser.parse_args=lambda self,*args,**kwargs: (_ for _ in ()).throw(Captured(self))
result={}
for name in sys.argv[1:]:
    try: runpy.run_path(str(Path(name).resolve()),run_name='__main__')
    except Captured as error:
        result[name]=[dict(flags=a.option_strings,dest=a.dest,default=a.default,choices=list(a.choices) if a.choices else None,required=a.required) for a in error.parser._actions if a.option_strings]
    except BaseException as error: result[name]=dict(error=str(error))
print(json.dumps(result,default=str))
