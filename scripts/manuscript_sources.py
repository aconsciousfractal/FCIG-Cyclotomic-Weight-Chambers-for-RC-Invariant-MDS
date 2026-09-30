"""Resolve the active TeX input tree; unused historical files are not manuscript."""
from pathlib import Path
import re


def active_tex_files(paper: Path) -> tuple[Path, ...]:
    ordered: list[Path] = []
    visiting: set[Path] = set()

    def visit(path: Path) -> None:
        path = path.resolve()
        if not path.is_relative_to(paper.resolve()):
            raise ValueError("tex_input_outside_paper")
        if path in visiting:
            raise ValueError("tex_input_cycle")
        visiting.add(path)
        ordered.append(path)
        source = re.sub(r"(?m)(?<!\\)%.*$", "", path.read_text(encoding="utf-8"))
        for name in re.findall(r"\\(?:input|include)\{([^}]+)\}", source):
            target = paper / name
            if not target.suffix:
                target = target.with_suffix(".tex")
            visit(target)
        visiting.remove(path)

    visit(paper / "main.tex")
    return tuple(ordered)
