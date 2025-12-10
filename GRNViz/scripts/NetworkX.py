import networkx as nx
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import random
from plotly.colors import qualitative
import sys
import pandas as pd
import traceback

def log_error(message, exception=None):
    """Print error message with optional exception details."""
    print(f"[ERROR] {message}", file=sys.stderr)
    if exception:
        print(f"[ERROR] Exception: {type(exception).__name__}: {str(exception)}", file=sys.stderr)
        print(f"[ERROR] Traceback:\n{traceback.format_exc()}", file=sys.stderr)

def log_info(message):
    """Print informational message."""
    print(f"[INFO] {message}")

def load_graph(file_path):
    G = nx.DiGraph()
    source_nodes = set()
    target_nodes = set()

    try:
        # Read data with automatic delimiter detection
        df = pd.read_csv(file_path, sep=None, engine="python", header=None)
        log_info(f"Successfully loaded file: {file_path} with {len(df)} rows")

        # Create NetworkX graph and nodes/edges
        G = nx.DiGraph()
        source_nodes = set()
        target_nodes = set()

        if len(df.columns) == 3:  # 3-column format (source, weight, target)
            df.columns = ["source", "weight", "target"]
            for _, row in df.iterrows():
                try:
                    weight = round(float(row["weight"]), 4)
                    G.add_edge(row["source"], row["target"], weight=weight)
                    source_nodes.add(row["source"])
                    target_nodes.add(row["target"])
                except (ValueError, TypeError) as e:
                    log_error(f"Invalid weight value at row: {row.to_dict()}", e)
                    continue
            log_info(f"Detected 3-column format. Created graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

        elif len(df.columns) == 2:  # 2-column format (gene, weight)
            df.columns = ["gene", "weight"]
            for _, row in df.iterrows():
                try:
                    weight = round(float(row["weight"]), 4)
                    G.add_node(row["gene"], weight=weight)
                    source_nodes.add(row["gene"])
                except (ValueError, TypeError) as e:
                    log_error(f"Invalid weight value at row: {row.to_dict()}", e)
                    continue
            log_info(f"Detected 2-column format. Created graph with {G.number_of_nodes()} nodes")
        else:
            log_error(f"Input file must have 2 or 3 columns, got {len(df.columns)} columns")
            return None, None, None

    except FileNotFoundError as e:
        log_error(f"File not found: {file_path}", e)
        return None, None, None
    except pd.errors.EmptyDataError as e:
        log_error(f"Input file is empty: {file_path}", e)
        return None, None, None
    except ValueError as e:
        log_error(f"Value error while processing file: {file_path}", e)
        return None, None, None
    except Exception as e:
        log_error(f"Unexpected error while loading graph from: {file_path}", e)
        return None, None, None

    return G, source_nodes, target_nodes

def sample_graph_by_degree(G, top_n=10):
    """Sample subgraph with top N nodes by degree."""
    try:
        degree_dict = dict(G.degree())
        if not degree_dict:  # Empty graph check
            log_info("Graph is empty, returning original graph")
            return G
        # Use the smaller value between actual graph size and requested top_n
        actual_n = min(top_n, len(degree_dict))
        if actual_n < top_n:
            log_info(f"Requested top {top_n} nodes, but only {actual_n} available")
        top_nodes = sorted(degree_dict, key=degree_dict.get, reverse=True)[:actual_n]
        subgraph = G.subgraph(top_nodes).copy()
        log_info(f"Sampled subgraph with {subgraph.number_of_nodes()} nodes")
        return subgraph
    except Exception as e:
        log_error("Failed to sample graph by degree", e)
        return G

def calculate_node_size(G, method="degree", scale=20, min_size=10):
    """Calculate node sizes based on centrality method."""
    try:
        if len(G.nodes()) == 0:  # Empty graph check
            log_info("Graph is empty, returning empty node sizes")
            return []
        if method == "degree":
            centrality = dict(G.degree())
        elif method == "betweenness":
            centrality = nx.betweenness_centrality(G)
        elif method == "pagerank":
            try:
                centrality = nx.pagerank(G)
            except nx.PowerIterationFailedConvergence as e:
                log_error("PageRank failed to converge, falling back to degree centrality", e)
                centrality = dict(G.degree())
        else:
            log_error(f"Unknown centrality method: {method}. Using 'degree' as fallback.")
            centrality = dict(G.degree())

        max_value = max(centrality.values()) or 1  # Avoid division by zero
        min_value = min(centrality.values())

        # Normalize node sizes between min_size and scale
        if max_value == min_value:
            # All nodes have the same centrality value
            return [min_size + (scale - min_size) / 2 for _ in centrality.values()]

        return [
            ((value - min_value) / (max_value - min_value)) * (scale - min_size) + min_size
            for value in centrality.values()
        ]
    except Exception as e:
        log_error("Failed to calculate node sizes", e)
        return [min_size for _ in G.nodes()]

def plot_graph_plotly(G, source_nodes, target_nodes, node_sizes):
    """Create Plotly visualization of the graph."""
    try:
        pos = nx.spring_layout(G, seed=42)

        edge_colors = random.choices(qualitative.Plotly, k=max(len(G.edges), 1))
        node_colors = random.choices(qualitative.Set3, k=max(len(G.nodes), 1))

        edge_trace = []
        for i, edge in enumerate(G.edges(data=True)):
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            weight = edge[2].get('weight', 1)
            edge_trace.append(go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                line=dict(width=weight, color=edge_colors[i]),
                hoverinfo='none',
                mode='lines'
            ))

        node_x = []
        node_y = []
        node_text = []
        for node in G.nodes:
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            text=node_text,
            mode='markers+text',
            textposition='top center',
            hoverinfo='text',
            marker=dict(
                showscale=False,
                colorscale='Viridis',
                color=node_colors,
                size=node_sizes,
                line_width=1
            )
        )

        fig = go.Figure(data=edge_trace + [node_trace],
                        layout=go.Layout(
                            title=dict(text='Gene Regulatory Network', font=dict(size=16)),
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            annotations=[dict(
                                text="",
                                showarrow=False,
                                xref="paper", yref="paper"
                            )],
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                        ))

        log_info("Successfully created Plotly visualization")
        return fig
    except Exception as e:
        log_error("Failed to create Plotly visualization", e)
        raise

def main():
    # Validate command line arguments
    if len(sys.argv) < 4:
        log_error(f"Insufficient arguments. Expected at least 3 arguments (input_file, output_file, top_n), got {len(sys.argv) - 1}")
        log_error("Usage: python NetworkX.py <input_file> <output_file> <top_n> [centrality_method]")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    # Validate top_n argument
    try:
        top_n = int(sys.argv[3])
        if top_n <= 0:
            log_error(f"top_n must be a positive integer, got {top_n}")
            sys.exit(1)
    except ValueError as e:
        log_error(f"Invalid top_n value: '{sys.argv[3]}'. Must be an integer.", e)
        sys.exit(1)

    centrality_method = sys.argv[4] if len(sys.argv) > 4 else "degree"
    valid_methods = ["degree", "betweenness", "pagerank"]
    if centrality_method not in valid_methods:
        log_error(f"Invalid centrality method: '{centrality_method}'. Valid options: {valid_methods}")
        log_info(f"Using default method: 'degree'")
        centrality_method = "degree"

    log_info(f"Processing input file: {input_file_path}")
    log_info(f"Output file: {output_file_path}")
    log_info(f"Top N nodes: {top_n}")
    log_info(f"Centrality method: {centrality_method}")

    G, source_nodes, target_nodes = load_graph(input_file_path)

    if G is None:
        log_error("Failed to load graph. Exiting.")
        sys.exit(1)

    if G.number_of_nodes() == 0:
        log_error("Graph has no nodes. Cannot create visualization.")
        sys.exit(1)

    try:
        sampled_G = sample_graph_by_degree(G, top_n=top_n)
        node_sizes = calculate_node_size(sampled_G, method=centrality_method, scale=80, min_size=10)

        fig = plot_graph_plotly(sampled_G, source_nodes, target_nodes, node_sizes)
        fig_json = fig.to_json()

        with open(output_file_path, "w") as json_file:
            json_file.write(fig_json)
        log_info(f"Successfully saved output to: {output_file_path}")

    except PermissionError as e:
        log_error(f"Permission denied when writing to: {output_file_path}", e)
        sys.exit(1)
    except Exception as e:
        log_error("Failed to create and save visualization", e)
        sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log_error("Unexpected error occurred", e)
        sys.exit(1)
