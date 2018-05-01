import json


def get_leaves(tree):
    if len(tree['children']) == 0:
        return 1
    leaves = 0
    for child in tree['children']:
        leaves += get_leaves(child)
    return leaves


def get_height(tree):
    if len(tree['children']) == 0:
        return 0
    max_height = 0
    for child in tree['children']:
        max_height = max(max_height, get_height(child))

    # For single-node unmutated clones
    if ((len(tree['children']) == 1 and
            len(tree['children'][0]['data']['mutations']) == 0)):
        return max_height
    return 1 + max_height


def get_height_muts(tree):
    max_height = 0
    for child in tree['children']:
        max_height = max(max_height, get_height_muts(child))
    return len(tree['data']['mutations']) + max_height


def get_node_num(tree):
    ttl = 0
    for child in tree['children']:
        ttl += get_node_num(child)
    return ttl + 1


def get_max_leaf_dist(tree):
    if len(tree['children']) == 0:
        return 0
    if len(tree['children']) == 1:
        return get_max_leaf_dist(tree['children'][0])
    widths = []
    heights = []
    for child in tree['children']:
        heights.append(get_height(child))
        widths.append(get_max_leaf_dist(child))
    heights = sorted(heights, reverse=True)
    return max(max(widths), heights[0] + heights[1] + 2)


def tree_compare(found, correct, error):
    found = json.loads(found)['tree']
    correct = json.loads(correct)['tree']
    metrics = [get_leaves, get_height, get_height_muts, get_node_num,
               get_max_leaf_dist]

    for metric_func in metrics:
        found_v = metric_func(found)
        correct_v = metric_func(correct)
        assert found_v == correct_v, \
            'Trees mismatch with {}: {} should be {}'.format(
                    metric_func.__name__, found_v, correct_v)
