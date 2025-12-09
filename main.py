from fastmcp import FastMCP
from typing import List, Optional, Dict, Tuple
import random

# ==============================================================================
# 1. Thermodynamic Parameters (SantaLucia J Jr. 1997) & Constants
# ==============================================================================
NN_PARAMS = {
    'aa': -1.0, 'tt': -1.0, 
    'at': -0.88, 
    'ta': -0.58, 
    'ca': -1.45, 'tg': -1.45, 
    'gt': -1.44, 'ac': -1.44, 
    'ct': -1.28, 'ag': -1.28, 
    'ga': -1.30, 'tc': -1.30, 
    'cg': -2.17, 
    'gc': -2.24, 
    'gg': -1.84, 'cc': -1.84, 
    'terGCp': 0.98, 
    'terATp': 1.03, 
    'penaltySelf': 0.43 
}

GLOBAL_BULGE = 3
GLOBAL_HAIRPIN_LOOP = 7
STEP_OFFSET = 2

# ==============================================================================
# 2. Helper Functions (String Manipulation & Bio-utils)
# ==============================================================================

def get_complement(seq: str) -> str:
    """Returns the complementary DNA sequence."""
    trans_table = str.maketrans("ATGC", "TACG")
    return seq.translate(trans_table)

def get_reverse_complement(seq: str) -> str:
    """Returns the reverse complementary DNA sequence."""
    return get_complement(seq)[::-1]

def calculate_gc(seq: str) -> float:
    """Calculates GC content percentage."""
    if not seq: return 0.0
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100

def get_slice(main_seq: str, parallel_seq: str, step: int) -> str:
    """
    Port of String.prototype.slicer from JS.
    Simulates sliding window alignment.
    """
    len_main = len(main_seq)
    len_para = len(parallel_seq)
    
    # Determine short/long and case c
    if len_main < len_para:
        short, long, c = len_main, len_para, 1
    elif len_main > len_para:
        short, long, c = len_para, len_main, 2
    else:
        short = long = len_main
        c = 3
    
    # Logic mapping based on JS implementation
    sliced_seq = ""
    
    if c == 1:
        if step <= short:
            sliced_seq = main_seq[0:step]
        elif short < step <= long:
            sliced_seq = main_seq[0:short]
        else: # step > long
            start = step - long
            end = start + (short - (step - long))
            sliced_seq = main_seq[start:end]
            
    elif c == 2:
        if step <= short:
            sliced_seq = main_seq[0:step]
        elif short < step <= long:
            start = step - short
            sliced_seq = main_seq[start:start+short]
        else:
            start = step - short
            length = short - (step - long)
            sliced_seq = main_seq[start:start+length]
            
    elif c == 3:
        if step <= short:
            sliced_seq = main_seq[0:step]
        elif short < step < 2 * short:
            start = step - short
            length = short - (step - short)
            sliced_seq = main_seq[start:start+length]
            
    return sliced_seq

def remove_segment(seq: str, index: int, length: int) -> str:
    """Port of String.prototype.insertStr for deletion (creating bulges)."""
    return seq[:index] + seq[index+length:]

# ==============================================================================
# 3. Core Algorithm: dG Calculation (The Heart of exCGen)
# ==============================================================================

def calculate_dg(seq_new: str, seq_old: str) -> float:
    """
    Calculates Gibbs Free Energy (dG) using NN model.
    seq_new: 5'->3', seq_old: 5'->3' (comparison target)
    """
    if not seq_new or not seq_old:
        return 0.0

    seq_new = seq_new.strip()
    seq_old = seq_old.strip()
    pairing_len = len(seq_new)
    
    # Convert comparison sequence from 3' to 5' effectively for pairing check
    # In JS: spOld_rev[pairingLen - (i + 1)] = spOld[i]
    # This reverses seq_old
    seq_old_rev = seq_old[::-1] 
    
    # Pairing logic
    # AT->a, TA->b, GC->c, CG->d, Others->x
    pairing_code = []
    
    for i in range(pairing_len):
        if i >= len(seq_old_rev): break # Safety check
        
        pair = seq_new[i] + seq_old_rev[i]
        if pair == "AT": code = 'a'
        elif pair == "TA": code = 'b'
        elif pair == "GC": code = 'c'
        elif pair == "CG": code = 'd'
        else: code = 'x'
        pairing_code.append(code)
        
    pair_result = "_" + "".join(pairing_code) + "_"
    
    dG = 0.0
    
    # Sliding window of size 2 over the code string
    for i in range(len(pair_result) - 1):
        nn_pair = pair_result[i:i+2]
        
        # Look up dG params
        # Handling the switch case from JS
        val = 0.0
        if nn_pair in ['aa', 'bb']: val = NN_PARAMS['aa'] # aa/tt are same (-1)
        elif nn_pair == 'ab': val = NN_PARAMS['at']
        elif nn_pair == 'ba': val = NN_PARAMS['ta']
        elif nn_pair in ['ad', 'da']: val = NN_PARAMS['ca'] # ca/tg
        elif nn_pair in ['cb', 'bc']: val = NN_PARAMS['gt'] # gt/ac
        elif nn_pair in ['db', 'bd']: val = NN_PARAMS['ct'] # ct/ag
        elif nn_pair in ['ca', 'ac']: val = NN_PARAMS['ga'] # ga/tc
        elif nn_pair == 'dc': val = NN_PARAMS['cg']
        elif nn_pair == 'cd': val = NN_PARAMS['gc']
        elif nn_pair in ['cc', 'dd']: val = NN_PARAMS['gg'] # gg/cc
        
        # Terminal penalties
        elif nn_pair in ['_c', '_d', 'c_', 'd_']: val = NN_PARAMS['terGCp']
        elif nn_pair in ['_a', '_b', 'a_', 'b_']: val = NN_PARAMS['terATp']
        
        dG += val
        
    return dG

# ==============================================================================
# 4. Generators & Gating Logic
# ==============================================================================

def generate_random_sequence(preset_pattern: str, gc_ratio: float, tolerance: float = 5.0) -> Optional[dict]:
    """Generates a random sequence matching preset and GC requirements."""
    # Base pools
    arr_ATGC = ["A", "T", "G", "C"]
    arr_ATATGC = ["A", "T", "A", "T", "A", "T", "G", "C", "G", "C"] # Low GC
    arr_ATGCGC = ["A", "T", "A", "T", "G", "C", "G", "C", "G", "C"] # High GC
    
    if gc_ratio <= 45: base_pool = arr_ATATGC
    elif gc_ratio > 55: base_pool = arr_ATGCGC
    else: base_pool = arr_ATGC
    
    max_attempts = 1000
    for _ in range(max_attempts):
        p_strand = []
        for char in preset_pattern:
            if char == '*':
                p_strand.append(random.choice(base_pool))
            elif char in ['A', 'T', 'G', 'C']:
                p_strand.append(char)
            # Handle R, Y, S, W if needed (omitted for brevity, can add if requested)
            else:
                p_strand.append(random.choice(base_pool))
        
        p_strand_str = "".join(p_strand)
        curr_gc = calculate_gc(p_strand_str)
        
        if gc_ratio - tolerance <= curr_gc <= gc_ratio + tolerance:
            # Valid GC
            minus_strand = get_complement(p_strand_str[::-1])
            duplex_dg = calculate_dg(p_strand_str, minus_strand)
            return {
                "plus": p_strand_str,
                "minus": minus_strand,
                "gcr": curr_gc,
                "duplex": duplex_dg
            }
            
    return None

def generate_bulges(seq: str, bulge_len: int = 3, start_pos: int = 3) -> List[str]:
    """Generates variations of sequence with missing bases (bulges)."""
    bulges = [seq] # Index 0 is original
    for h in range(1, bulge_len + 1):
        # JS loop: i <= inputSeq.length - (h + 2 * startPos)
        limit = len(seq) - (h + 2 * start_pos)
        for i in range(limit + 1):
            bulged = remove_segment(seq, i + start_pos, h)
            bulges.append(bulged)
    return bulges

def check_hairpin(plus_strand: str, threshold: float) -> bool:
    """Checks for hairpin formation."""
    # Simplified loop based on JS logic
    start_pos = 3
    limit_h = GLOBAL_HAIRPIN_LOOP
    
    lowest_dg = 0.0
    
    # Loop for hairpin stem length
    for h in range(3, limit_h + 1):
        limit_i = len(plus_strand) - (h + 2 * start_pos)
        for i in range(limit_i + 1):
            # JS: slice="on" logic (simplified)
            left_part = plus_strand[0 : i + start_pos]
            right_part = plus_strand[i + start_pos + h :]
            
            # Slicer usage in JS hairpinMaker is complex. 
            # Logic: Compare Left vs Right for self-annealing
            # We use the slicer logic implicitly by checking interaction
            
            # Note: JS implementation creates specific left/right objects using slicer.
            # Here we approximate the heavy lifting of `cal_dG` on self parts.
            # Using the slicer on these parts:
            
            # Re-implementing specific JS hairpin loop logic:
            # It compares right part against left part
            
            # For strict port, we need the slicer helper again here? 
            # Actually, `cal_dG` expects 5'-3' strings.
            # Let's trust `cal_dG` handles length mismatch via internal logic if passed directly?
            # No, `cal_dG` iterates `pairingLen`.
            
            # Let's use the explicit logic:
            sliced_right = get_slice(right_part, left_part, len(left_part)) # Step = short length
            sliced_left = get_slice(left_part, right_part, len(right_part))
            
            dg = calculate_dg(sliced_right, sliced_left) + NN_PARAMS['penaltySelf']
            if dg < lowest_dg:
                lowest_dg = dg
                
    if lowest_dg < threshold: # More stable than threshold (negative values)
        return False # Reject
    return True

def check_cross_dimer(candidate: dict, existing_list: List[dict], 
                      threshold_cross: float, threshold_selected: float,
                      selected_list: List[dict]) -> bool:
    """
    Checks interaction against:
    1. Self (Self-dimer)
    2. Previously generated sequences (Cross-dimer)
    3. Selected (Input) sequences (Selected-dimer)
    """
    
    new_plus = candidate['plus']
    new_minus = candidate['minus']
    
    # 1. Self Dimer Check
    # Generate bulges for self
    bulges_plus = generate_bulges(new_plus, GLOBAL_BULGE)
    bulges_minus = generate_bulges(new_minus, GLOBAL_BULGE)
    
    lowest_self = 0.0
    
    # Iterate bulges
    for b_seq in bulges_plus:
        # Sliding window (JS: j = StepOffset ... )
        limit = len(b_seq) - STEP_OFFSET
        for j in range(STEP_OFFSET, limit):
            # new_plus vs itself (with bulge)
            dg = calculate_dg(get_slice(b_seq, new_plus, j), get_slice(new_plus, b_seq, j)) + NN_PARAMS['penaltySelf']
            if dg < lowest_self: lowest_self = dg
            
    if lowest_self < threshold_cross: return False
    
    # 2. Cross Dimer (vs Existing Generated)
    lowest_cross = 0.0
    for exist in existing_list:
        exist_plus = exist['plus']
        exist_minus = exist['minus']
        
        # Check Candidate Bulges vs Existing
        for b_seq in bulges_plus:
             limit = len(b_seq) - STEP_OFFSET
             for j in range(STEP_OFFSET, limit):
                 # vs Plus
                 dg_pp = calculate_dg(get_slice(b_seq, exist_plus, j), get_slice(exist_plus, b_seq, j))
                 if dg_pp < lowest_cross: lowest_cross = dg_pp
                 # vs Minus
                 dg_pm = calculate_dg(get_slice(b_seq, exist_minus, j), get_slice(exist_minus, b_seq, j))
                 if dg_pm < lowest_cross: lowest_cross = dg_pm
                 
        if lowest_cross < threshold_cross: return False

    # 3. Selected Dimer (vs Input constraints)
    lowest_selected = 0.0
    for sel in selected_list:
        sel_seq = sel['seq']
        # Check Candidate Bulges vs Selected Sequence
        for b_seq in bulges_plus:
            # Step length logic from JS
            step_len = max(len(b_seq), len(sel_seq))
            limit = step_len - STEP_OFFSET
            
            for j in range(STEP_OFFSET, limit):
                 dg = calculate_dg(get_slice(b_seq, sel_seq, j), get_slice(sel_seq, b_seq, j))
                 if dg < lowest_selected: lowest_selected = dg
        
        # Adjust threshold for length difference (JS Logic)
        local_threshold = threshold_selected
        if len(sel_seq) < len(new_plus):
            local_threshold = (threshold_selected * len(sel_seq) / len(new_plus))
            
        if lowest_selected < local_threshold: return False

    return True

# ==============================================================================
# 5. MCP Server Definition
# ==============================================================================

mcp = FastMCP("exCGen-DNA-Architect")

@mcp.tool()
def design_orthogonal_dna(
    length: int,
    count: int,
    gc_min: float = 40.0,
    gc_max: float = 60.0,
    preset_pattern: str = None,
    existing_sequences: list[str] = None,
    avoid_sequences: list[str] = None
) -> str:
    """
    Generates a set of orthogonal (mutually exclusive) DNA sequences using the SantaLucia NN thermodynamic model.
    
    Args:
        length: Length of DNA (bp).
        count: Number of pairs to generate.
        gc_min: Minimum GC content (%).
        gc_max: Maximum GC content (%).
        preset_pattern: Pattern string (e.g., "****ATGC****"). '*' is random. Length must match 'length'.
        existing_sequences: List of DNA strings that new sequences must NOT bind to (for Request #2).
        avoid_sequences: List of DNA strings to strictly avoid (e.g., restriction sites).
    """
    
    # Input validation / Setup
    if preset_pattern and len(preset_pattern) != length:
        return f"Error: Preset pattern length ({len(preset_pattern)}) does not match requested length ({length})."
    
    if not preset_pattern:
        preset_pattern = "*" * length

    # Parse inputs
    selected_objs = []
    if existing_sequences:
        for s in existing_sequences:
            s = s.upper().strip()
            # Calculate simplistic dG self for obj structure
            selected_objs.append({'seq': s, 'id': 'input'})

    avoid_list = avoid_sequences if avoid_sequences else []
    
    # Initialize Storage
    passed_sequences = []
    attempts = 0
    max_total_attempts = count * 1000 # Safety break
    
    results_log = []
    
    # Main Loop
    while len(passed_sequences) < count and attempts < max_total_attempts:
        attempts += 1
        
        # 1. Distribute Target GC
        # Linearly distribute GC between min and max based on current progress
        if count > 1:
            step = (gc_max - gc_min) / (count - 1)
            target_gc = gc_min + (step * len(passed_sequences))
        else:
            target_gc = (gc_min + gc_max) / 2
            
        # 2. Generate Candidate
        candidate = generate_random_sequence(preset_pattern, target_gc)
        if not candidate:
            continue
            
        # 3. Calculate Dynamic Thresholds
        # JS: threshold = (percent * duplex_dG) / 100
        # Default percentages from JS
        th_hair_per = 10
        th_cross_per = 25
        th_sel_per = 30
        
        duplex_dg = candidate['duplex']
        th_hair = (th_hair_per * duplex_dg) / 100
        th_cross = (th_cross_per * duplex_dg) / 100
        th_sel = (th_sel_per * duplex_dg) / 100
        
        # 4. Gating (Filters)
        
        # A. Explicit Avoids (String matching)
        rejected = False
        for avoid in avoid_list:
            if avoid in candidate['plus'] or avoid in candidate['minus']:
                rejected = True
                break
        if rejected: continue

        # B. Hairpin Check
        if not check_hairpin(candidate['plus'], th_hair):
            continue
            
        # C. Cross/Self/Selected Check
        if not check_cross_dimer(candidate, passed_sequences, th_cross, th_sel, selected_objs):
            continue
            
        # 5. Success
        passed_sequences.append(candidate)
        results_log.append(f"Seq {len(passed_sequences)}: {candidate['plus']} (GC: {candidate['gcr']:.1f}%, dG: {candidate['duplex']:.2f})")

    # Format Output
    output_text = f"Successfully generated {len(passed_sequences)}/{count} orthogonal pairs.\n"
    output_text += "-" * 40 + "\n"
    
    for idx, obj in enumerate(passed_sequences):
        output_text += f"Pair {idx+1}:\n"
        output_text += f"  Plus : {obj['plus']}\n"
        output_text += f"  Minus: {obj['minus']}\n"
        output_text += f"  GC: {obj['gcr']:.1f}%, Duplex dG: {obj['duplex']:.2f} kcal/mol\n\n"
        
    return output_text

# ==============================================================================
# Run Server
# ==============================================================================
if __name__ == "__main__":
    mcp.run()