def get_bonded_pairs(sequance, structure):
    pairs = []
    basis =[]

    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")":
            j = stack.pop()
            pairs.append((j, i))
            basis.append((sequence[i], sequence[j]))
            
    return pairs, basis