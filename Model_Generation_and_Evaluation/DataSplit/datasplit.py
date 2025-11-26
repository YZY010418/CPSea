import argparse
import csv
import random
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Cluster-based train/validation split for a TSV file.')
    parser.add_argument('-i', '--input_tsv', type=str, required=True,
                        help='Input TSV file with cluster_center and member columns.')
    parser.add_argument('--valid_ratio', type=float, default=0.02,
                        help='Proportion of unique clusters to sample for the validation set (default: 0.02).')
    args = parser.parse_args()

    # (1) Read the TSV file and group members by cluster_center
    center_to_members = defaultdict(list)
    with open(args.input_tsv, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        center_col_index = 0
        member_col_index = 1
        
        for row in reader:
            if len(row) > 1:
                center = row[center_col_index]
                member = row[member_col_index]
                center_to_members[center].append(member)

    # (2) Randomly sample cluster centers
    all_centers = list(center_to_members.keys())
    total_clusters = len(all_centers)
    num_valid_clusters = int(total_clusters * args.valid_ratio)
    valid_centers_sampled = set(random.sample(all_centers, num_valid_clusters))

    # (3) Split members and save to train/valid files
    train_members = []
    valid_members = []
    
    train_clusters_final = []
    valid_clusters_final = []

    for center, members in center_to_members.items():
        if center in valid_centers_sampled:
            valid_members.extend(members)
            valid_clusters_final.append(center)
        else:
            train_members.extend(members)
            train_clusters_final.append(center)

    with open('CPCore_valid.txt', 'w') as f:
        f.write('\n'.join(valid_members) + '\n')

    with open('CPCore_train.txt', 'w') as f:
        f.write('\n'.join(train_members) + '\n')

    # (4) Print summary statistics
    actual_valid_clusters = len(valid_clusters_final)
    actual_train_clusters = len(train_clusters_final)
    total_members = len(train_members) + len(valid_members)
    
    actual_cluster_ratio = actual_valid_clusters / total_clusters if total_clusters > 0 else 0
    actual_member_ratio = len(valid_members) / total_members if total_members > 0 else 0

    print("\n--- Sampling Summary ---")
    print(f"Total Unique Clusters: {total_clusters}")
    print(f"Requested Validation Ratio: {args.valid_ratio:.4f}")
    
    print("\n[Clusters]")
    print(f"  Train Clusters: {actual_train_clusters}")
    print(f"  Valid Clusters: {actual_valid_clusters}")
    print(f"  Actual Valid Cluster Ratio: {actual_cluster_ratio:.4f}")
    
    print("\n[Members]")
    print(f"  Total Members: {total_members}")
    print(f"  Valid Members: {len(valid_members)}")
    print(f"  Actual Valid Member Ratio: {actual_member_ratio:.4f}")

if __name__ == '__main__':
    main()
