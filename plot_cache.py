import os


def reset_cache(dir):
    # print("reset")
    for file in os.listdir(dir):
        path = os.path.join(dir, file)

        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            reset_cache(path)
            os.rmdir(path)


if __name__ == "__main__":
    reset_cache("plot_cache")
