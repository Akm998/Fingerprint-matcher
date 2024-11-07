import math
from collections import defaultdict

def euclidean_distance(minutiae1, minutiae2):
    """
    Calculate the Euclidean distance between two minutiae points.
    """
    x1, y1, _ = minutiae1
    x2, y2, _ = minutiae2
    
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def angular_difference(minutiae1, minutiae2):
    """
    Calculate the angular difference between two minutiae points.
    """
    _, _, theta1 = minutiae1
    _, _, theta2 = minutiae2
    
    angle_diff = abs(theta1 - theta2) % 360
    if angle_diff > 180:
        angle_diff = 360 - angle_diff
        
    return angle_diff

def load_xyt_file(filepath):
    """
    Load minutiae data from a .xyt file.
    """
    minutiae_points = []
    with open(filepath, 'r') as file:
        for line in file:
            values = list(map(int, line.split()))
            minutiae_points.append((values[0], values[1], values[2]))
    return minutiae_points

def create_pairwise_table(minutiae_points):
    """
    Create a pairwise minutiae comparison table (local structures).
    """
    pairwise_table = {}
    
    for i, minutiae1 in enumerate(minutiae_points):
        for j, minutiae2 in enumerate(minutiae_points):
            if i != j:
                distance = euclidean_distance(minutiae1, minutiae2)
                angle_diff = angular_difference(minutiae1, minutiae2)
                pairwise_table[(i + 1, j + 1)] = (distance, angle_diff)
    
    return pairwise_table

def match_pairwise_tables(table1, table2, dist_threshold=20, angle_threshold=15):
    """
    Match pairwise minutiae tables from two fingerprints using the Bozorth3 algorithm.
    """
    matched_pairs = defaultdict(list)
    
    for pair1, (dist1, angle1) in table1.items():
        if pair1 in table2:
            dist2, angle2 = table2[pair1]
            
            if abs(dist1 - dist2) <= dist_threshold and abs(angle1 - angle2) <= angle_threshold:
                matched_pairs[pair1].append((dist1, angle1))
    
    print(f"Total matched pairs: {len(matched_pairs)}")  # Debugging print
    return matched_pairs

def cluster_matching_pairs(matched_pairs, dist_threshold=20, angle_threshold=15):
    """
    Cluster matched minutiae pairs into connected components (graph-based clustering).
    """
    clusters = []
    visited = set()
    adjacency_list = defaultdict(list)

    # Build adjacency list: add edges between pairs that share a minutiae point and are geometrically consistent
    for pair1 in matched_pairs:
        for pair2 in matched_pairs:
            if pair1 != pair2:
                # Check if they share a minutiae point AND are geometrically consistent
                if (pair1[0] == pair2[0] or pair1[1] == pair2[1]) and _is_geometrically_consistent(matched_pairs[pair1], matched_pairs[pair2], dist_threshold, angle_threshold):
                    adjacency_list[pair1].append(pair2)
                    adjacency_list[pair2].append(pair1)

    # Perform DFS to find connected components
    def dfs(pair, cluster):
        stack = [pair]
        while stack:
            current_pair = stack.pop()
            if current_pair not in visited:
                visited.add(current_pair)
                cluster.add(current_pair)
                for neighbor in adjacency_list[current_pair]:
                    if neighbor not in visited:
                        stack.append(neighbor)

    # Traverse all matched pairs to create clusters
    for pair in matched_pairs:
        if pair not in visited:
            cluster = set()
            dfs(pair, cluster)
            clusters.append(cluster)
    
    print(f"Total clusters: {len(clusters)}")  # Debugging print
    return clusters

def _is_geometrically_consistent(pair1_data, pair2_data, dist_threshold, angle_threshold):
    """
    Check if two matched minutiae pairs are geometrically consistent.
    """
    (dist1, angle1) = pair1_data[0]
    (dist2, angle2) = pair2_data[0]

    return abs(dist1 - dist2) <= dist_threshold and abs(angle1 - angle2) <= angle_threshold

def bozorth3_score(file1, file2, dist_threshold=20, angle_threshold=15):
    """
    Calculate the Bozorth3 matching score between two fingerprint minutiae sets.
    """
    # Load minutiae points from both files
    minutiae_points1 = load_xyt_file(file1)
    minutiae_points2 = load_xyt_file(file2)
    
    # Create pairwise tables for both fingerprints
    pairwise_table1 = create_pairwise_table(minutiae_points1)
    pairwise_table2 = create_pairwise_table(minutiae_points2)
    
    # Match the pairwise tables to find matched pairs
    matched_pairs = match_pairwise_tables(pairwise_table1, pairwise_table2, dist_threshold, angle_threshold)
    
    # Cluster the matched pairs based on mutual consistency
    clusters = cluster_matching_pairs(matched_pairs, dist_threshold, angle_threshold)
    
    # The score is the number of clusters found
    return len(clusters)

# Example usage
file1 = "minutiae_file1.xyt"
file2 = "minutiae_file2.xyt"

score = bozorth3_score(file1, file2)
print(f"Bozorth3 Matching Score: {score}")
