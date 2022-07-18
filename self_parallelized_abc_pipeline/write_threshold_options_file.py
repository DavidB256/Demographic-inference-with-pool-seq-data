
threshold = "0.001"
repeats = 10

with open("threshold_options.txt", "w") as f:
    for i in range(repeats):
        f.write(threshold + "\n")
