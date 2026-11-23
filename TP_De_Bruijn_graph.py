from collections import defaultdict
from collections import Counter

filepath = input("What the file path of your file ? : ").strip().strip("'").strip('"')

kmer_size = int(input("What k-mer size do you wan to use ? : ").strip())


def read_fastq(path):
    reads = []
    with open(path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            reads.append(sequence)
            print("read read")
    return reads

def get_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1 ):
        kmers.append(sequence[i:i+k])
    return kmers

def build_de_bruijn_graph(reads, k):
    graph = defaultdict(list)
    nodes = set()

    for read in reads:
        for kmer in get_kmers(read, k):
            left_node = kmer[:-1]
            right_node = kmer[1:]
            graph[left_node].append(right_node)
            nodes.add(left_node)
            nodes.add(right_node)

    return graph, nodes

def find_start_node(graph, nodes):
    in_deg = Counter()
    out_deg = Counter()

    for u, neighbors in graph.items():
        out_deg[u] += len(neighbors)
        for v in neighbors:
            in_deg[v] += 1

    for node in nodes:
        if out_deg[node] - in_deg[node] == 1 :
            return node
    return next(iter(graph))

def eulerian_path(graph, start_node):
    local_graph = {u: list(v) for u,v in graph.items()}
    stack = [start_node]
    path = []

    while stack:
        u = stack[-1]
        if local_graph[u]:
            v = local_graph[u].pop()
            stack.append(v)
        else:
            path.append(stack.pop())
    path.reverse()
    return path

def export_to_dot(graph, output_file = "graph.dot"):
    with open(output_file, 'w') as f:
        f.write("digraph G {\n")
        f.write("    rankdir=LR;\n")
        f.write("    node [shape=box, style=rounded, fontname=\"Courier\"];\n")
        for u, neighbors in graph.items():
            for v in set(neighbors):
                weight = neighbors.count(v)
                f.write(f'"{u}" -> "{v}" [label="{weight}x"];\n')
        f.write("}\n")
    print(f"Graph exported successfully : {output_file}")

def assemble_sequence(graph, start_node):
    path = eulerian_path(graph, start_node)
    return path[0] + "".join(node[-1] for node in path[1:])

if __name__ == "__main__":
    reads = read_fastq(filepath)
    print(f"{len(reads)} séquences lues.")

    graph, nodes = build_de_bruijn_graph(reads, kmer_size)
    print(f"Graphe construit : {len(nodes)} nœuds.")

    output_dot = r"TP_De_Bruijn_graph_project\Data\debruijn_graph.dot"
    export_to_dot(graph, output_dot)

    simple_graph = {u: list(set(v)) for u, v in graph.items()}
    start = find_start_node(simple_graph, nodes)
    contig = assemble_sequence(simple_graph, start)

    print(f"\nNœud de départ : {start}")
    print(f"Séquence assemblée ({len(contig)} pb) :\n{contig}")


