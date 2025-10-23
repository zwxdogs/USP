from probe.single_probe import Single_probe

s_p = Single_probe(5e6, 0.6, [0, 0, 0], 0.01, 0.015)

line = s_p.get_line()

print("Line equations:")
for i, l in enumerate(line):
    print(f"Line {i+1}: {l[0]}*x + {l[1]}*y + {l[2]} = 0")
