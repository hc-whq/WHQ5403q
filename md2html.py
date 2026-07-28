import sys
import re

def convert_md_to_html(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    html = ["<html><body>"]
    in_list = False

    for line in lines:
        line = line.strip()
        if not line:
            if in_list:
                html.append("</ul>")
                in_list = False
            continue

        # Headers
        if line.startswith("# "):
            html.append(f"<h1>{line[2:]}</h1>")
        elif line.startswith("## "):
            html.append(f"<h2>{line[3:]}</h2>")
        elif line.startswith("### "):
            html.append(f"<h3>{line[4:]}</h3>")
        # List items
        elif line.startswith("- "):
            if not in_list:
                html.append("<ul>")
                in_list = True
            content = line[2:]
            # Bold
            content = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', content)
            html.append(f"<li>{content}</li>")
        else:
            if in_list:
                html.append("</ul>")
                in_list = False
            # Regular paragraph
            content = line
            # Bold
            content = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', content)
            html.append(f"<p>{content}</p>")

    if in_list:
        html.append("</ul>")
    
    html.append("</body></html>")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(html))

if __name__ == "__main__":
    convert_md_to_html(sys.argv[1], sys.argv[2])
