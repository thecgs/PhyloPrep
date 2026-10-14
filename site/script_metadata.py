"""Collect CLI-declared versions and arguments without executing script work."""
import argparse
import json
from pathlib import Path
import runpy
import sys


class ParserCaptured(BaseException):
    def __init__(self, parser):
        self.parser = parser


def collect_script_metadata(paths):
    result = {}
    original = argparse.ArgumentParser.parse_args
    def capture(parser, *_args, **_kwargs):
        raise ParserCaptured(parser)
    argparse.ArgumentParser.parse_args = capture
    try:
        for item in paths:
            path = Path(item).resolve()
            sys.path.insert(0, str(path.parent))
            try:
                namespace = runpy.run_path(str(path), run_name='__main__')
                result[path.name] = dict(version=namespace.get('__version__'), actions=[])
            except ParserCaptured as captured:
                actions = captured.parser._actions
                version = next((action.version for action in actions if isinstance(action, argparse._VersionAction)), None)
                result[path.name] = dict(version=version, actions=[
                    dict(flags=action.option_strings, dest=action.dest, default=action.default,
                         choices=list(action.choices) if action.choices else None, required=action.required)
                    for action in actions if action.option_strings
                ])
            finally:
                sys.path.pop(0)
    finally:
        argparse.ArgumentParser.parse_args = original
    return result


if __name__ == '__main__':
    print(json.dumps(collect_script_metadata(sys.argv[1:]), default=str))
