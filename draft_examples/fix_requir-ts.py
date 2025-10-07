input_file = "requirements.txt"          # your corrupted one
output_file = "requirements_clean.txt"   # new clean version

with open(input_file, "rb") as f:
    data = f.read()

# Decode as UTF-8 ignoring errors, then remove null bytes and control chars
text = data.decode("utf-8", errors="ignore")
text = text.replace("\x00", "").replace("\r", "")

# Keep only printable characters and newlines
text = "".join(ch for ch in text if ch.isprintable() or ch in "\n\t")

# Write clean UTF-8 text
with open(output_file, "w", encoding="utf-8", newline="\n") as f:
    f.write(text)

print(f"✅ Clean UTF-8 file written to {output_file}")
