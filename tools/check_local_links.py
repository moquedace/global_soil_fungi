"""Check that relative links in the Markdown files point to existing files."""
import pathlib
import re
import sys

root = pathlib.Path(__file__).resolve().parent.parent
link = re.compile(r"\[[^\]]*\]\(([^)#\s]+)(?:#[^)]*)?\)")
missing = []
for md in root.rglob("*.md"):
    if ".git" in md.parts:
        continue
    for target in link.findall(md.read_text(encoding="utf-8")):
        if re.match(r"^(https?|mailto):", target):
            continue
        if not (md.parent / target).resolve().exists():
            missing.append(f"{md.relative_to(root)} -> {target}")
if missing:
    print("Broken local links:\n" + "\n".join(missing))
    sys.exit(1)
print("All local links resolve.")
