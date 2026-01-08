import pandas as pd
import numpy as np
import itertools, operator


def GenerateTMCircleRelativeCenters(df, membraneThickness, R):
    start_center = (0, 0)
    TMCircles_min = df[df["IMO"]=="TMhelix"]["length"].min()
    membraneThickness_max = (TMCircles_min-1)*2*R
    if membraneThickness > membraneThickness_max:
        raise ValueError(f"membraneThickness must be below {membraneThickness_max}.")
    
    TMUnits_idx_centers = {}
    for idx, row in df.iterrows():
        centers = []
        if row["IMO"] != "TMhelix":
            continue

        N_circles = row["length"]
        L = 2*R*N_circles
        sinTheta = membraneThickness / L
        cosTheta = np.sqrt(1.0 - sinTheta**2)
        centers += [start_center]
        for i in range(N_circles-1): # next N-1 circles
            x, y = centers[-1]
            sign = 1 if (i % 2 == 0) else -1
            dx = sign * 2.0 * R * cosTheta 
            dy = 2.0 * R * sinTheta     

            next_center = (x + dx, y + dy)
            centers.append(next_center)
        TMUnits_idx_centers[idx] = centers

    return TMUnits_idx_centers
        


def Convert_dff32df(deeptmhmm_gff3_fp):
    df = pd.read_csv(deeptmhmm_gff3_fp, sep="\t", comment="#", header=None,
                names=["name", "IMO", "start", "end"], 
                usecols=[0,1,2,3])
    df["length"] = df["end"] - df["start"] + 1
    return df

def Convert_3line2df(deeptmhmm_3line_fp):
    # get header/seq/IMO
    with open(deeptmhmm_3line_fp) as f:
        lines = [line.strip() for line in f if line.strip()]
    header = lines[0]         
    seq = lines[1]           
    imo = lines[2] 
    assert len(seq) == len(imo)
    # position,amino acid(AA),IMO columns
    n = len(seq)
    df = pd.DataFrame({
        "position": range(1, n + 1),
        "AA": list(seq),
        "IMO": list(imo)
    })

    df["segment_id"] = (df["IMO"] != df["IMO"].shift()).cumsum()
    seg = df.groupby("segment_id").agg(
        start=("segment_id", lambda x: x.index.min() + 1),   
        end=("segment_id", lambda x: x.index.max() + 1),
        length=("segment_id", "size")
        )

    df = df.merge(seg, on="segment_id").drop(columns="segment_id")
    return df

def MoveCoords(coords, new_start):
    """
    Move a set of (x,y) to new positions.
    The first or last (x,y) is moved to new_start while keeping all other positions fixed relative to each other.
    """
    X, Y = new_start      # new desired position of first point
    x0, y0 = coords[0]

    dx = X - x0
    dy = Y - y0
    return [(x + dx, y + dy) for x, y in coords]

def _pick_k_with_replacement_targetSUM_sorted(arr, k, target_sum=None):
    """
    Pick k items (with replacement) from 1D array, and the sum of picked K items is target_sum
    """
    n_items = len(arr)
    combinations = []
    for idx_combo in itertools.combinations_with_replacement(range(n_items), k):
        items = [arr[idx] for idx in idx_combo]
        if target_sum is not None and sum(items) != target_sum:
            continue
        item_min = min(items)
        item_max = max(items)
        var = item_max - item_min

        combinations.append({"var": var,
                        "items": items,
                        "indices": idx_combo,
                        })

    combinations.sort(key=lambda item: (item["var"], item["items"]))
    return combinations

def CirclesCombinations(N, max_a, max_b):
    if N <= 3: # N_min = 3
        return N
    else:
        x_max = (N-1)//2
        y_max = (N+1)//4
        
        a_range = np.arange(1,max_a+1)
        
        y_range = np.arange(1,y_max+1)
        b_range = np.arange(1,max_b+1)
        b2_range = 2*b_range

        
        for y in y_range:
            b2_item_pos_pool = _pick_k_with_replacement_targetSUM_sorted(b2_range, k=y, target_sum=None)
            for b2_item_pos_dict in b2_item_pos_pool:
                b2_values = b2_item_pos_dict["items"]
                b2_values_sum = np.sum(b2_values)
                N_rest = N - b2_values_sum
                x = 2*y-1
                a_val_pos_pool_dict = _pick_k_with_replacement_targetSUM_sorted(a_range, k=x, target_sum=N_rest)
                if len(a_val_pos_pool_dict) == 0:
                    continue
                else:
                    a_values = a_val_pos_pool_dict[0]["items"]
                    bSet = [b2/2 for b2 in b2_values for _ in range(2)]
                    ans_list = [v for pair in zip(bSet, a_values) for v in pair] + [bSet[-1]]
                    return ans_list
                

def _next_circle_center(prev_center, R, upward, degree):
    x, y = prev_center
    if upward:
        angle_deg=degree
    else:
        angle_deg=-degree
    angle_rad = np.deg2rad(angle_deg)
    dx = 2 * R * np.cos(angle_rad)
    dy = 2 * R * np.sin(angle_rad)
    return (x + dx, y + dy)

def _gen_straight_nCircleCenters(pre_center, n, R, upward):
    centers = [pre_center]
    for _ in range(int(n)):
        pre_center = centers[-1]
        c_next = _next_circle_center(pre_center, R, upward, degree=90)
        centers.append(c_next)
    return centers 

def _gen_peak_nCircleCenters(pre_center, n, R, upward):
    if n == 1:
        center = _next_circle_center(pre_center, R, upward, degree=30)
        bridge_center = _next_circle_center(center, R, upward=operator.not_(upward), degree=30)
        return [center, bridge_center]
    elif n ==2:
        center1 = _next_circle_center(pre_center, R, upward, degree=60)
        center2 = _next_circle_center(center1, R, upward, degree=0)
        bridge_center = _next_circle_center(center2, R, upward=operator.not_(upward), degree=60)
        return [center1, center2, bridge_center]



def _gen_horizontalS(pre_center, R, restCircles_comb_is_list, upward):
    centers_list = [pre_center]
    for idx, val in enumerate(restCircles_comb_is_list):
        if idx%4 in (2,3):
            upward = not upward

        if idx%4 == 0:
            pre_center = centers_list[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center, val-1, R, upward)
            centers_list += next_CircleCenters[1:]
        elif idx%4 == 1:
            pre_center = centers_list[-1]
            peakcenters_bridge = _gen_peak_nCircleCenters(pre_center, val, R, upward)
            centers_list += peakcenters_bridge
        elif idx%4 == 2:
            pre_center = centers_list[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center, val-1, R, upward)
            centers_list += next_CircleCenters[1:]
        else: # idx%4 == 3
            pre_center = centers_list[-1]
            peakcenters_bridge = _gen_peak_nCircleCenters(pre_center, val, R, upward)
            centers_list += peakcenters_bridge
    return centers_list


def _peak_curve_centers123(pre_center, R, upward, remainingCircles_comb):
    # peak curve
    peak123_bridge_centers = []
    if remainingCircles_comb == 1:
        center = center = _next_circle_center(pre_center, R, upward, degree=30)
        bridge_center = _next_circle_center(center, R, upward=operator.not_(upward), degree=30)
        peak123_bridge_centers += [center, bridge_center]
    elif remainingCircles_comb == 2:
        # remainingCircles_comb == 2:
        center1 = _next_circle_center(pre_center, R, upward, degree=60)
        center2 = _next_circle_center(center1, R, upward, degree=0)
        bridge_center = _next_circle_center(center2, R, upward=operator.not_(upward), degree=60)
        peak123_bridge_centers += [center1, center2, bridge_center]
    else:
        # remainingCircles_comb == 3:
        center1 = _next_circle_center(pre_center, R, upward, degree=90)
        center2 = _next_circle_center(center1, R, upward, degree=30)
        center3 = _next_circle_center(center2, R, operator.not_(upward), degree=60)
        bridge_center = _next_circle_center(center3, R, upward=operator.not_(upward), degree=90)
        peak123_bridge_centers += [center1, center2, center3, bridge_center]
    return peak123_bridge_centers




def _over0_on_check_remainingCircles(pre_center, length, away, R, upward, max_a, max_b):
    centers_bridge = []
    remainingCircles = length - away * 2
    remainingCircles_comb = CirclesCombinations(remainingCircles, max_a, max_b)
    if isinstance(remainingCircles_comb, list):
        # left away side
        left_start_center = pre_center
        left_next_CircleCenters = _gen_straight_nCircleCenters(left_start_center, away+1, R, upward)
        # peak curve
        horizontalSCenters = _gen_horizontalS(left_next_CircleCenters[-1], R, remainingCircles_comb, upward)
        # right away side
        right_start_center = horizontalSCenters[-1]
        right_next_CircleCenters = _gen_straight_nCircleCenters(right_start_center, away+1, R, operator.not_(upward))
        centers_bridge = centers_bridge + left_next_CircleCenters + horizontalSCenters + right_next_CircleCenters
        return centers_bridge
    else:
        # remaining circles is limited in [1,2,3]
        # left away side
        left_start_center = pre_center
        left_next_CircleCenters = _gen_straight_nCircleCenters(left_start_center, away, R, upward)
        # peak curve
        peak123_curve_centers = _peak_curve_centers123(pre_center=left_next_CircleCenters[-1], R=R, remainingCircles_comb=remainingCircles_comb, upward=upward)
        # right
        right_start_center = peak123_curve_centers[-1]
        right_next_CircleCenters = _gen_straight_nCircleCenters(right_start_center, away, R, operator.not_(upward))
        centers_bridge = centers_bridge + left_next_CircleCenters + peak123_curve_centers + right_next_CircleCenters
        return centers_bridge



def _equal0_on_check_remainingCircles(pre_center, away, R, upward):
    centers_bridge = []
    left_next_CircleCenters = []
    left_start_center = pre_center
    for _ in range(away):
        left_next_CircleCenters.append(_next_circle_center(left_start_center, R, upward, degree=30))

    right_start_center = _next_circle_center(left_next_CircleCenters[-1], R, upward, degree=0) 
    right_next_CircleCenters = []
    for _ in range(away-1):
        right_next_CircleCenters.append(_next_circle_center(right_start_center, R, upward=operator.not_(upward), degree=30))
    bridge_center = _next_circle_center(right_next_CircleCenters[-1], R, upward=operator.not_(upward), degree=90)
    centers_bridge = centers_bridge + left_next_CircleCenters + right_next_CircleCenters + [bridge_center]
    return centers_bridge



def _below0_on_check_remainingCircles(pre_center, length, R, upward, max_a, max_b):
    centers_bridge = []
    remainingCircles_comb = CirclesCombinations(length, max_a, max_b)
    if isinstance(remainingCircles_comb, list):
        # left away side
        left_start_center = pre_center
        left_next_CircleCenters = _gen_straight_nCircleCenters(left_start_center, 1, R, upward)
        # peak curve
        horizontalSCenters = _gen_horizontalS(left_next_CircleCenters[-1], R, remainingCircles_comb, upward)
        # right away side
        right_start_center = horizontalSCenters[-1]
        right_next_CircleCenters = _gen_straight_nCircleCenters(right_start_center, 1, R, operator.not_(upward))
        centers_bridge = centers_bridge + left_next_CircleCenters + horizontalSCenters + right_next_CircleCenters
        return centers_bridge
    else:
        # remaining circles is limited in [1,2,3]
        # peak curve
        peak123_curve_centers = _peak_curve_centers123(pre_center=pre_center, R=R, remainingCircles_comb=remainingCircles_comb, upward=upward)
        centers_bridge = centers_bridge + peak123_curve_centers
        return centers_bridge




def _not_extracellular_not_Nterm_not_Cterm(pre_center, length, away, R, max_a, max_b):
    upward = False
    check_remainingCircles = length - away * 2
    if check_remainingCircles > 0 :
        return _over0_on_check_remainingCircles(pre_center, length, away, R, upward, max_a, max_b)
    elif check_remainingCircles < 0 :
        return _below0_on_check_remainingCircles(pre_center, length, R, upward, max_a, max_b)
    else:
        return _equal0_on_check_remainingCircles(pre_center, away, R, upward)
    

def _extracellular_not_Nterm_not_Cterm(pre_center, length, away, R, max_a, max_b):
    upward = True
    check_remainingCircles = length - away * 2
    if check_remainingCircles > 0 :
        return _over0_on_check_remainingCircles(pre_center, length, away, R, upward, max_a, max_b)
    elif check_remainingCircles < 0 :
        return _below0_on_check_remainingCircles(pre_center, length, R, upward, max_a, max_b)
    else:
        return _equal0_on_check_remainingCircles(pre_center, away, R, upward)



def AddExtracellularNotTMCenters(pre_center, length, away, R, max_a, max_b):
    return _extracellular_not_Nterm_not_Cterm(pre_center, length, away, R, max_a, max_b)

def AddIntracellularNotTMCenters(pre_center, length, away, R, max_a, max_b):
    return _not_extracellular_not_Nterm_not_Cterm(pre_center, length, away, R, max_a, max_b)


###################################
# Add Nterm
###################################
def AddNterm_Centers(pre_center, length, away, R, IMO, max_a, max_b):
    if IMO == "inside":
        return _add_Intracellular_Nterm_Centers(pre_center, length, away, R, max_a, max_b)
    else:
        return _add_Extracellular_Nterm_Centers(pre_center, length, away, R, max_a, max_b)

def _add_Intracellular_Nterm_Centers(pre_center, length, away, R, max_a, max_b):
    # Nterm inside
    centers_bridge = [pre_center]
    remainingCircles = length - away # one side

    if length <= away:
        pre_center = centers_bridge[-1]
        next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=True)
        centers_bridge += next_CircleCenters[1:]
        return centers_bridge
    else:
        restCircles_comb = CirclesCombinations(remainingCircles, max_a, max_b)
        if isinstance(restCircles_comb, list):
            Scurve_centers = _gen_horizontalS(pre_center=centers_bridge[-1], R=R, restCircles_comb_is_list=restCircles_comb, upward=False)
            centers_bridge += Scurve_centers[1:]
            # right side away
            pre_center = centers_bridge[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=away, R=R, upward=True)
            centers_bridge += next_CircleCenters[1:]
            # bridge center
            bridge_center = _next_circle_center(centers_bridge[-1], R, upward=True, degree=90)
            centers_bridge += [bridge_center]
            return centers_bridge
        else:
            # short Nterm
            pre_center = centers_bridge[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=True)
            centers_bridge += next_CircleCenters[1:]
            return centers_bridge
    

def _add_Extracellular_Nterm_Centers(pre_center, length, away, R, max_a, max_b):
    # Nterm inside
    centers_bridge = [pre_center]
    remainingCircles = length - away # one side
    if length <= away:
        pre_center = centers_bridge[-1]
        next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=False)
        centers_bridge += next_CircleCenters[1:]
        return centers_bridge
    else:
        restCircles_comb = CirclesCombinations(remainingCircles, max_a, max_b)
        if isinstance(restCircles_comb, list):
            Scurve_centers = _gen_horizontalS(pre_center=centers_bridge[-1], R=R, restCircles_comb_is_list=restCircles_comb, upward=True)
            centers_bridge += Scurve_centers[1:]
            # right side away
            pre_center = centers_bridge[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=away, R=R, upward=False)
            centers_bridge += next_CircleCenters[1:]
            # bridge center
            bridge_center = _next_circle_center(centers_bridge[-1], R, upward=False, degree=90)
            centers_bridge += [bridge_center]
            return centers_bridge
        else:
            # short Nterm
            pre_center = centers_bridge[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=False)
            centers_bridge += next_CircleCenters[1:]
            return centers_bridge

    


###################################
# Add Cterm
###################################
def AddCterm_Centers(pre_center, length, away, R, IMO, max_a, max_b):
    if IMO == "inside":
        return _add_Intracellular_Cterm_Centers(pre_center, length, away, R, max_a, max_b)
    else:
        return _add_Extracellular_Cterm_Centers(pre_center, length, away, R, max_a, max_b)
    

def _add_Intracellular_Cterm_Centers(pre_center, length, away, R, max_a, max_b):
    # Cterm inside
    centers_bridge = [pre_center]
    remainingCircles = length - away # one side

    if length <= away:
        pre_center = centers_bridge[-1]
        next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=False)
        centers_bridge += next_CircleCenters[1:]
        return centers_bridge
    else:
        restCircles_comb = CirclesCombinations(remainingCircles, max_a, max_b)
        if isinstance(restCircles_comb, list):
            # left one side
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=away+1, R=R, upward=False)
            # S curve
            Scurve_centers = _gen_horizontalS(pre_center=next_CircleCenters[-1], R=R, restCircles_comb_is_list=restCircles_comb, upward=False)
            centers_bridge = centers_bridge + next_CircleCenters[1:] + Scurve_centers[1:]
            return centers_bridge
        else:
            # short Nterm
            pre_center = centers_bridge[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=False)
            centers_bridge += next_CircleCenters[1:]
            return centers_bridge
        

def _add_Extracellular_Cterm_Centers(pre_center, length, away, R, max_a, max_b):
    # Cterm outside
    centers_bridge = [pre_center]
    remainingCircles = length - away # one side

    if length <= away:
        pre_center = centers_bridge[-1]
        next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=True)
        centers_bridge += next_CircleCenters[1:]
        return centers_bridge
    else:
        restCircles_comb = CirclesCombinations(remainingCircles, max_a, max_b)
        if isinstance(restCircles_comb, list):
            # left one side
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=away+1, R=R, upward=True)
            # S curve
            Scurve_centers = _gen_horizontalS(pre_center=next_CircleCenters[-1], R=R, restCircles_comb_is_list=restCircles_comb, upward=True)
            centers_bridge = centers_bridge + next_CircleCenters[1:] + Scurve_centers[1:]
            return centers_bridge
        else:
            # short Nterm
            pre_center = centers_bridge[-1]
            next_CircleCenters = _gen_straight_nCircleCenters(pre_center=pre_center, n=length, R=R, upward=True)
            centers_bridge += next_CircleCenters[1:]
            return centers_bridge