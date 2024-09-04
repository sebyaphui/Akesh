def print_dict_with_function(d, f):
    for k, v in d.items():
        print(f"{f(k)} : {v}")
