"""
Restructure gene map coordinates and groupings.

File: restructure_gene_map.py
Author: Akaash Venkat, Audi Liu
Updated: 2026 (Python 3 modernization)

This module reorganizes gene positions and groupings in the SVG visualization
using circular layout algorithms and group classification heuristics.
"""

import logging
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration
INFO_DIR = Path("info_files")
INPUT_DIR = Path("input_files")
SVG_DIR = Path("svg_files")

GENE_DATABASE_FILE = INFO_DIR / "gene_database.txt"
INTERMEDIATE_GENES_FILE = INFO_DIR / "intermediate_genes.txt"
GENE_COORDS_FILE = INFO_DIR / "gene_coords.txt"
GENE_GROUP_FILE = INFO_DIR / "gene_group.txt"
GENE_GROUPINGS_FILE = INFO_DIR / "gene_groupings.txt"
GROUPINGS_FILE = INPUT_DIR / "grouping_details.txt"
BASE_SVG = SVG_DIR / "original_gene_map.svg"
RESTRUCTURED_SVG = SVG_DIR / "restructured_gene_map.svg"

# Global data structures
GENE_LIST: List[Tuple[str, Dict[str, float]]] = []
GENE_GROUP: Dict[str, List[str]] = {}
GENE_GROUPINGS: Dict[int, List[str]] = {}
B_D_PAIR: Dict[str, str] = {}


def read_database() -> None:
    """Parse gene database file."""
    if not GENE_DATABASE_FILE.exists():
        return

    try:
        with open(GENE_DATABASE_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split(' - ', 1)
                if len(parts) != 2:
                    continue

                gene_name = parts[0].strip()
                neighbors = {}

                for neighbor_entry in parts[1].split(', '):
                    match = re.match(r'(\w+)\(([\d.]+)\)', neighbor_entry.strip())
                    if match:
                        neighbor_name, confidence = match.groups()
                        neighbors[neighbor_name] = float(confidence)

                GENE_LIST.append((gene_name, neighbors))

        logger.info(f"Loaded {len(GENE_LIST)} genes from database")
    except Exception as e:
        logger.error(f"Error reading database: {e}")


def read_gene_groups() -> None:
    """Parse gene groups file."""
    if not GENE_GROUP_FILE.exists():
        return

    try:
        with open(GENE_GROUP_FILE, 'r') as f:
            current_group = None
            for line in f:
                line = line.strip()
                if not line or '---' in line:
                    continue

                if 'Group' in line and ':' in line:
                    current_group = line.split('Group ')[1].split(':')[0]
                    GENE_GROUP[current_group] = []
                elif current_group:
                    GENE_GROUP[current_group].append(line)

        logger.info(f"Loaded {len(GENE_GROUP)} gene groups")
    except Exception as e:
        logger.error(f"Error reading gene groups: {e}")


def read_gene_groupings() -> None:
    """Parse gene groupings file."""
    if not GROUPINGS_FILE.exists():
        return

    try:
        with open(GROUPINGS_FILE, 'r') as f:
            for line in f:
                if 'Group' in line and ':' in line:
                    grouping_id = int(line.split('Group ')[1].split(':')[0])
                    grouping_elements = line.split(': ', 1)[1].split(', ')
                    GENE_GROUPINGS[grouping_id] = grouping_elements

        logger.info(f"Loaded {len(GENE_GROUPINGS)} groupings")
    except Exception as e:
        logger.error(f"Error reading groupings: {e}")


def read_intermediate_pairs() -> None:
    """Parse intermediate gene pairs."""
    if not INTERMEDIATE_GENES_FILE.exists():
        return

    try:
        with open(INTERMEDIATE_GENES_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if line and ' : ' in line and 'following' not in line.lower():
                    b_gene, d_gene = line.split(' : ')
                    B_D_PAIR[b_gene] = d_gene

        logger.info(f"Loaded {len(B_D_PAIR)} B-D pairs")
    except Exception as e:
        logger.error(f"Error reading intermediate pairs: {e}")


def write_gene_coords(new_pos_dict: Dict[str, Tuple[str, str]]) -> None:
    """Write gene coordinates to file."""
    try:
        with open(GENE_COORDS_FILE, 'w') as f:
            f.write("Gene coordinates in restructured SVG file.\n\n")
            for gene, coords in sorted(new_pos_dict.items()):
                f.write(f"{gene} : [{coords[0]}, {coords[1]}]\n")
        logger.info(f"Wrote coordinates for {len(new_pos_dict)} genes")
    except Exception as e:
        logger.error(f"Error writing gene coords: {e}")


def write_gene_groupings(groupings: List[List[str]]) -> None:
    """Write reorganized gene groupings."""
    try:
        with open(GENE_GROUPINGS_FILE, 'w') as f:
            for i, group_genes in enumerate(groupings, 1):
                f.write(f"Group {i}:\n---\n")
                for gene in group_genes:
                    f.write(f"{gene}\n")
                f.write("\n\n\n")
        logger.info("Gene groupings file written")
    except Exception as e:
        logger.error(f"Error writing gene groupings: {e}")


def classify_genes(
    group1: List[str], group2: List[str],
    group3: List[str], group4: List[str]
) -> List[List[str]]:
    """Recursively classify unclassified genes based on connections."""
    groups = [group1[:], group2[:], group3[:], group4[:]]
    a_genes = set(GENE_GROUP.get("A", []))

    # Remove already classified genes
    for group in groups:
        a_genes -= set(group)

    # Iteratively classify remaining genes
    iteration = 0
    while a_genes and iteration < 100:
        iteration += 1
        classified = set()

        for gene in list(a_genes):
            # Find gene in GENE_LIST
            gene_neighbors = {}
            for name, neighbors in GENE_LIST:
                if name == gene:
                    gene_neighbors = neighbors
                    break

            if not gene_neighbors:
                continue

            # Count connections to each group
            counts = [0, 0, 0, 0]
            confidences = [0.0, 0.0, 0.0, 0.0]

            for neighbor, confidence in gene_neighbors.items():
                for i, group in enumerate(groups):
                    if neighbor in group:
                        counts[i] += 1
                        confidences[i] += confidence
                        break

            # Classify based on max count
            max_count = max(counts) if sum(counts) > 0 else 0
            if max_count == 0:
                continue

            max_indices = [i for i, c in enumerate(counts) if c == max_count]

            if len(max_indices) == 1:
                groups[max_indices[0]].append(gene)
                classified.add(gene)
            else:
                # Tie-breaker: use highest confidence
                best_idx = max(
                    max_indices,
                    key=lambda i: confidences[i]
                )
                groups[best_idx].append(gene)
                classified.add(gene)

        a_genes -= classified
        if not classified:
            break  # No progress made

    logger.info(f"Classification complete after {iteration} iterations")
    return groups


def store_positions(
    new_pos_dict: Dict[str, Tuple[str, str]],
    text_pos_dict: Dict[str, Tuple[str, str]],
    gene_list: List[str],
    a_d_map: Dict[str, List[str]],
    d_b_map: Dict[str, str],
    center_x: int,
    center_y: int,
    offset1: float,
    offset2: float,
    factor: float,
    radius_factor: int,
) -> None:
    """Calculate and store circular positions for genes."""
    radius = len(gene_list) * radius_factor

    for i, gene in enumerate(gene_list):
        angle = 2 * math.pi * i / len(gene_list)
        x = int(round(center_x + radius * math.cos(angle)))
        y = int(round(center_y + radius * math.sin(angle)))

        new_pos_dict[gene] = (str(x), str(y))
        text_pos_dict[gene] = (
            str(int(offset1 * x - offset2 * center_x)),
            str(int(offset1 * y - offset2 * center_y)),
        )

        # Position connecting D genes
        if gene in a_d_map:
            offset1 *= factor
            offset2 = offset1 - 1

            d_genes = a_d_map[gene]
            if len(d_genes) == 1:
                _position_d_gene(
                    d_genes[0], d_b_map, new_pos_dict, text_pos_dict,
                    center_x, center_y, radius, i, len(gene_list), offset1, offset2
                )
            else:
                for k, d_gene in enumerate(d_genes):
                    sub_angle = i - 0.25 + k * 0.5 / (len(d_genes) - 1)
                    _position_d_gene(
                        d_gene, d_b_map, new_pos_dict, text_pos_dict,
                        center_x, center_y, radius, sub_angle, len(gene_list), offset1, offset2
                    )

            offset1 /= factor
            offset2 = offset1 - 1


def _position_d_gene(
    d_gene: str,
    d_b_map: Dict[str, str],
    new_pos_dict: Dict[str, Tuple[str, str]],
    text_pos_dict: Dict[str, Tuple[str, str]],
    center_x: int,
    center_y: int,
    radius: int,
    angle_idx: float,
    total_genes: int,
    offset1: float,
    offset2: float,
) -> None:
    """Helper to position a D gene and its associated B gene."""
    radius_d = radius * 1.2
    angle_rad = 2 * math.pi * angle_idx / total_genes
    d_x = int(round(center_x + radius_d * math.cos(angle_rad)))
    d_y = int(round(center_y + radius_d * math.sin(angle_rad)))

    new_pos_dict[d_gene] = (str(d_x), str(d_y))
    text_pos_dict[d_gene] = (
        str(int(offset1 * d_x - offset2 * center_x)),
        str(int(offset1 * d_y - offset2 * center_y)),
    )

    if d_gene in d_b_map:
        b_gene = d_b_map[d_gene]
        radius_b = radius_d * 1.1
        b_x = int(round(center_x + radius_b * math.cos(angle_rad)))
        b_y = int(round(center_y + radius_b * math.sin(angle_rad)))

        new_pos_dict[b_gene] = (str(b_x), str(b_y))
        text_pos_dict[b_gene] = (
            str(int(offset1 * b_x - offset2 * center_x)),
            str(int(offset1 * b_y - offset2 * center_y)),
        )


def update_dict(d_gene: str, dictionary: Dict, gene_info: Tuple) -> None:
    """Update gene count and confidence tracking."""
    dictionary[d_gene][0] += 1
    if gene_info[1].get(d_gene, 0) > dictionary[d_gene][2]:
        dictionary[d_gene][2] = gene_info[1][d_gene]
        dictionary[d_gene][1] = gene_info[0]


def modify_svg_content(
    content: List[str],
    old_pos_dict: Dict[str, str],
    new_pos_dict: Dict[str, Tuple[str, str]],
    text_pos_dict: Dict[str, Tuple[str, str]],
    *groups: List[str],
) -> List[str]:
    """Update SVG with new positions and styling."""
    modified_content = []

    for i, line in enumerate(content):
        # Remove old styling
        if "font-size:" in line and "px;" in line:
            continue
        if 'r="20"' in line:
            line = line.replace('r="20"', 'r="4"')

        # Update edges with high confidence
        if 'line class="nw_edge"' in line:
            if '.0" stroke=' in line:
                continue  # Skip zero-confidence edges

            # Update coordinates for edges with non-zero confidence
            if ('.1" stroke=' in line or '.2" stroke=' in line):
                line = _update_edge_coordinates(line, old_pos_dict, new_pos_dict)

        # Update circle nodes
        if '<circle cx' in line:
            line = _update_circle(line, old_pos_dict, new_pos_dict)
        elif '<circle class' in line:
            line = _update_circle_class(line, old_pos_dict, new_pos_dict)

        # Update text labels
        if '<text ' in line:
            line = _update_text(line, old_pos_dict, new_pos_dict, text_pos_dict, groups)

        modified_content.append(line)

    return modified_content


def _update_edge_coordinates(line: str, old_pos: Dict, new_pos: Dict) -> str:
    """Update edge coordinates in SVG line elements."""
    try:
        # Extract old coordinates
        matches = re.findall(r'(x1|y1|x2|y2)="([\d.]+)"', line)
        if len(matches) < 4:
            return line

        coords = {m[0]: float(m[1]) for m in matches}
        old_key = f"{int(coords['x1'] - 0.5)} {int(coords['y1'] - 0.5)}"

        if old_key not in old_pos:
            return line

        gene1_name = old_pos[old_key]
        old_key2 = f"{int(coords['x2'] - 0.5)} {int(coords['y2'] - 0.5)}"

        if old_key2 not in old_pos:
            return line

        gene2_name = old_pos[old_key2]

        if gene1_name in new_pos and gene2_name in new_pos:
            # Update stroke properties
            match_width = re.search(r'stroke-width="([\d.]+)"', line)
            match_opacity = re.search(r'stroke-opacity="([\d.]+)"', line)

            if match_width and match_opacity:
                new_width = 0.45 * float(match_width.group(1))
                new_opacity = 0.45 * float(match_opacity.group(1))

                new_x1 = float(new_pos[gene1_name][0]) + 0.5
                new_y1 = float(new_pos[gene1_name][1]) + 0.5
                new_x2 = float(new_pos[gene2_name][0]) + 0.5
                new_y2 = float(new_pos[gene2_name][1]) + 0.5

                line = re.sub(
                    r'stroke-opacity="[\d.]+"',
                    f'stroke-opacity="{new_opacity}"',
                    line
                )
                line = re.sub(
                    r'stroke-width="[\d.]+"',
                    f'stroke-width="{new_width}"',
                    line
                )
                line = re.sub(r'x1="[\d.]+"', f'x1="{new_x1}"', line)
                line = re.sub(r'y1="[\d.]+"', f'y1="{new_y1}"', line)
                line = re.sub(r'x2="[\d.]+"', f'x2="{new_x2}"', line)
                line = re.sub(r'y2="[\d.]+"', f'y2="{new_y2}"', line)
    except Exception as e:
        logger.debug(f"Error updating edge coordinates: {e}")

    return line


def _update_circle(line: str, old_pos: Dict, new_pos: Dict) -> str:
    """Update circle node coordinates."""
    try:
        match_x = re.search(r'cx="([\d.]+)"', line)
        match_y = re.search(r'cy="([\d.]+)"', line)

        if match_x and match_y:
            old_key = f"{match_x.group(1)} {match_y.group(1)}"
            if old_key in old_pos:
                gene = old_pos[old_key]
                if gene in new_pos:
                    line = re.sub(r'cx="[\d.]+"', f'cx="{new_pos[gene][0]}"', line)
                    line = re.sub(r'cy="[\d.]+"', f'cy="{new_pos[gene][1]}"', line)
    except Exception as e:
        logger.debug(f"Error updating circle: {e}")

    return line


def _update_circle_class(line: str, old_pos: Dict, new_pos: Dict) -> str:
    """Update circle class coordinates (same as regular circle)."""
    return _update_circle(line, old_pos, new_pos)


def _update_text(
    line: str,
    old_pos: Dict,
    new_pos: Dict,
    text_pos: Dict,
    groups: Tuple[List[str]]
) -> str:
    """Update text element position and font size."""
    try:
        match_x = re.search(r'x="([\d.]+)"', line)
        match_y = re.search(r'y="([\d.]+)"', line)

        if match_x and match_y:
            old_text_x = float(match_x.group(1))
            old_text_y = float(match_y.group(1))
            old_key = f"{int(old_text_x - 18)} {int(old_text_y + 18)}"

            if old_key in old_pos:
                gene = old_pos[old_key]
                if gene in new_pos and gene in text_pos:
                    # Determine font size based on group
                    font_size = 17
                    for i, group in enumerate(groups):
                        if gene in group:
                            font_size = [20, 15, 20, 18][i]
                            break

                    new_text_x = text_pos[gene][0]
                    new_text_y = text_pos[gene][1]

                    line = re.sub(r'text-anchor="[^"]+"', 'text-anchor="middle"', line)
                    line = re.sub(r'x="[\d.]+"', f'x="{new_text_x}"', line)
                    line = re.sub(r'y="[\d.]+"', f'y="{new_text_y}"', line)
                    line = re.sub(r'font-size="[\d]+"', f'font-size="{font_size}"', line)

    except Exception as e:
        logger.debug(f"Error updating text: {e}")

    return line


def modify_base_svg(groups_expanded: List[List[str]]) -> None:
    """Load base SVG and apply structural modifications."""
    try:
        with open(BASE_SVG, 'r') as f:
            content = f.readlines()

        # Update SVG header
        content[0] = '<svg class="notselectable" height="5000" id="svg_network_image" width="3500" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n'

        old_pos_dict = {}
        new_pos_dict = {}
        text_pos_dict = {}

        GROUP1, GROUP2, GROUP3, GROUP4 = groups_expanded
        d_list = GENE_GROUP.get("D", [])
        a_list = GENE_GROUP.get("A", [])

        # Build A-D mapping
        A_D_maps = [{}, {}, {}, {}]
        for d_gene in d_list:
            for i, d_map in enumerate(A_D_maps):
                d_map[d_gene] = [0, "", 0.0]

        for gene_name, neighbors in GENE_LIST:
            a_gene = gene_name
            for d_gene in d_list:
                if d_gene in neighbors:
                    if a_gene in GROUP1:
                        update_dict(d_gene, A_D_maps[0], (gene_name, neighbors))
                    elif a_gene in GROUP2:
                        update_dict(d_gene, A_D_maps[1], (gene_name, neighbors))
                    elif a_gene in GROUP3:
                        update_dict(d_gene, A_D_maps[2], (gene_name, neighbors))
                    elif a_gene in GROUP4:
                        update_dict(d_gene, A_D_maps[3], (gene_name, neighbors))

        # Determine best A gene for each D gene
        A_D_final = {a: [] for a in a_list}
        for d_gene in d_list:
            counts = [A_D_maps[i][d_gene][0] for i in range(4)]
            max_count = max(counts) if counts else 0
            max_indices = [i for i, c in enumerate(counts) if c == max_count]

            if len(max_indices) == 1:
                a_gene = A_D_maps[max_indices[0]][d_gene][1]
                A_D_final[a_gene].append(d_gene)
            else:
                confidences = [A_D_maps[i][d_gene][2] for i in max_indices]
                best_idx = max_indices[confidences.index(max(confidences))]
                a_gene = A_D_maps[best_idx][d_gene][1]
                A_D_final[a_gene].append(d_gene)

        # Create D-B mapping
        D_B_PAIR = {d: b for b, d in B_D_PAIR.items()}

        # Apply manual repositioning rules
        _apply_repositioning_rules(GROUP1, GROUP2)

        # Store positions for each group
        store_positions(new_pos_dict, text_pos_dict, GROUP1, A_D_final, D_B_PAIR, 1750, 400, 1.45, 0.45, 1.34, 12)
        store_positions(new_pos_dict, text_pos_dict, GROUP2, A_D_final, D_B_PAIR, 1750, 1600, 1.05, 0.05, 0.99, 11)
        store_positions(new_pos_dict, text_pos_dict, GROUP3, A_D_final, D_B_PAIR, 1750, 2900, 1.25, 0.25, 1.11, 12)
        store_positions(new_pos_dict, text_pos_dict, GROUP4, A_D_final, D_B_PAIR, 1750, 3950, 1.09, 0.09, 0.95, 12)

        # Extract old positions
        for i, line in enumerate(content):
            if 'g class="nwnodecontainer"' in line and 'data-safe_div_label' in line:
                match_gene = re.search(r'data-safe_div_label="([^"]+)"', line)
                match_x = re.search(r'data-x_pos="([\d.]+)"', line)
                match_y = re.search(r'data-y_pos="([\d.]+)"', line)

                if match_gene and match_x and match_y:
                    gene_name = match_gene.group(1)
                    x_pos = match_x.group(1)
                    y_pos = match_y.group(1)
                    old_pos_dict[f"{x_pos} {y_pos}"] = gene_name

        # Apply modifications
        content = modify_svg_content(content, old_pos_dict, new_pos_dict, text_pos_dict, GROUP1, GROUP2, GROUP3, GROUP4)

        # Write output
        with open(RESTRUCTURED_SVG, 'w') as f:
            f.writelines(content)

        write_gene_coords(new_pos_dict)
        logger.info("SVG restructuring completed")

    except Exception as e:
        logger.error(f"Error modifying SVG: {e}")


def _apply_repositioning_rules(group1: List[str], group2: List[str]) -> None:
    """Apply manual repositioning rules for specific genes."""
    try:
        # Reposition genes in Group 1
        for gene in ["RBP3", "FDFT1"]:
            if gene in group1:
                group1.remove(gene)
        group1.insert(2, "RBP3")
        group1.insert(4, "FDFT1")

        # Reposition genes in Group 2
        for gene in ["EMC1", "ALB"]:
            if gene in group2:
                group2.remove(gene)
        group2.append("EMC1")
        group2.insert(51, "ALB")

        reposition_genes = ["ALB", "ENO4", "NEK2", "AKT1", "RHO", "PDE6A", "GAPDH", "STAT3"]
        for gene in reposition_genes:
            if gene in group2:
                group2.remove(gene)
        for gene in reposition_genes:
            group2.insert(len(group2) - 1, gene)

        reposition_genes = ["SC5D", "MSMO1", "HSD17B7", "CYP51A1", "RCVRN"]
        for gene in reposition_genes:
            if gene in group2:
                group2.remove(gene)
        for gene in reposition_genes:
            group2.insert(len(group2) - 16, gene)

    except Exception as e:
        logger.warning(f"Error applying repositioning rules: {e}")


def main() -> None:
    """Main execution."""
    try:
        # Load data
        read_database()
        read_gene_groups()
        read_gene_groupings()
        read_intermediate_pairs()

        # Classify genes
        groups = classify_genes(
            GENE_GROUPINGS.get(1, []),
            GENE_GROUPINGS.get(2, []),
            GENE_GROUPINGS.get(3, []),
            GENE_GROUPINGS.get(4, [])
        )

        # Write groupings and modify SVG
        write_gene_groupings(groups)
        modify_base_svg(groups)

        logger.info("Gene map restructuring completed successfully")

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        raise


if __name__ == "__main__":
    main()
