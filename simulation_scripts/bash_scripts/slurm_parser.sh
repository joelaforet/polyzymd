#!/bin/bash

# SLURM Node Resource Parser - Fixed Version
# Processes 'scontrol show nodes' output for visualization and resource optimization

# Function to extract value from key=value pairs
extract_value() {
    local line="$1"
    local key="$2"
    # For AllocTRES and CfgTRES, we need to capture everything until the next key
    if [[ "$key" == "AllocTRES" || "$key" == "CfgTRES" ]]; then
        echo "$line" | sed -n "s/.*${key}=\([^[:space:]]*\).*/\1/p"
    else
        echo "$line" | grep -o "${key}=[^[:space:]]*" | cut -d'=' -f2-
    fi
}

# Fixed extract_gpu_info function
extract_gpu_info() {
    local gres="$1"
    local gpu_type="none"
    local gpu_total=0
    
    # Handle empty or missing gres
    if [[ -z "$gres" || "$gres" == "N/A" || "$gres" == "(null)" ]]; then
        echo "none,0"
        return
    fi
    
    # Handle multiple GPU types on same node (sum them up)
    if [[ "$gres" =~ gpu: ]]; then
        local temp_gres="$gres"
        local total=0
        local primary_type=""
        
        # Extract all gpu:type:count patterns
        while [[ "$temp_gres" =~ gpu:([^:,]+):([0-9]+) ]]; do
            local type="${BASH_REMATCH[1]}"
            local count="${BASH_REMATCH[2]}"
            
            # Use the first/primary GPU type for display
            if [[ -z "$primary_type" ]]; then
                primary_type="$type"
            fi
            
            total=$((total + count))
            temp_gres="${temp_gres#*${BASH_REMATCH[0]}}"
        done
        
        if [[ $total -gt 0 && -n "$primary_type" ]]; then
            gpu_type="$primary_type"
            gpu_total=$total
        fi
    fi
    
    # Ensure we always return valid values
    if [[ -z "$gpu_type" ]]; then
        gpu_type="none"
    fi
    if [[ -z "$gpu_total" ]]; then
        gpu_total=0
    fi
    
    echo "${gpu_type},${gpu_total}"
}

# Improved extract_alloc_gpus function
extract_alloc_gpus() {
    local alloc_tres="$1"
    local gpu_alloc=0
    
    # Debug: uncomment the next line to see the AllocTRES field
    # echo "DEBUG AllocTRES: $alloc_tres" >&2
    
    # Handle empty/missing AllocTRES
    if [[ -z "$alloc_tres" || "$alloc_tres" == "N/A" ]]; then
        echo "0"
        return
    fi
    
    # Look for various GPU allocation formats
    # Format: gres/gpu=number
    if [[ "$alloc_tres" =~ gres/gpu=([0-9]+) ]]; then
        gpu_alloc="${BASH_REMATCH[1]}"
    # Format: gres/gpu:type=number (handles complex types like h100_3g.40gb)
    elif [[ "$alloc_tres" =~ gres/gpu:[^=,]+=[0-9]+([^,]*,gres/gpu:[^=,]+=([0-9]+))* ]]; then
        # Sum all GPU allocations for complex cases
        local total=0
        local temp_tres="$alloc_tres"
        while [[ "$temp_tres" =~ gres/gpu:[^=,]+=([0-9]+) ]]; do
            total=$((total + ${BASH_REMATCH[1]}))
            temp_tres="${temp_tres#*${BASH_REMATCH[0]}}"
        done
        gpu_alloc=$total
    # Format: gres/gpu:type:number (alternative format)
    elif [[ "$alloc_tres" =~ gres/gpu:[^:,]+:([0-9]+) ]]; then
        gpu_alloc="${BASH_REMATCH[1]}"
    # Format: gres/gpu:number
    elif [[ "$alloc_tres" =~ gres/gpu:([0-9]+) ]]; then
        gpu_alloc="${BASH_REMATCH[1]}"
    fi
    
    echo "$gpu_alloc"
}

# Cache the scontrol output to avoid multiple calls
SCONTROL_OUTPUT=""
get_scontrol_output() {
    if [[ -z "$SCONTROL_OUTPUT" ]]; then
        SCONTROL_OUTPUT=$(scontrol show nodes)
    fi
    echo "$SCONTROL_OUTPUT"
}

# Main processing function
process_nodes() {
    local output_format="$1"
    
    # Header for CSV output
    if [[ "$output_format" == "csv" ]]; then
        echo "NodeName,State,CPUTotal,CPUAlloc,CPUFree,CPUUtil,MemTotal,MemAlloc,MemFree,MemUtil,GPUType,GPUTotal,GPUAlloc,GPUFree,GPUUtil,Partitions,Features"
    elif [[ "$output_format" == "json" ]]; then
        echo "["
        local first_node=true
    fi
    
    # Process nodes line by line, building complete node records
    local current_node=""
    local in_node=false
    
    while IFS= read -r line; do
        # Check if this line starts a new node
        if [[ "$line" =~ ^NodeName= ]]; then
            # Process the previous node if we have one
            if [[ -n "$current_node" ]]; then
                process_single_node "$current_node" "$output_format" "$first_node"
                first_node=false
            fi
            # Start new node
            current_node="$line"
            in_node=true
        elif [[ "$in_node" == true ]]; then
            # Continue building current node record
            current_node="$current_node $line"
        fi
    done < <(get_scontrol_output)
    
    # Process the last node
    if [[ -n "$current_node" ]]; then
        process_single_node "$current_node" "$output_format" "$first_node"
    fi
    
    if [[ "$output_format" == "json" ]]; then
        echo "]"
    fi
}

# Function to process a single node record
process_single_node() {
    local node_block="$1"
    local output_format="$2"
    local first_node="$3"
    
    # Extract basic node information
    local node_name=$(extract_value "$node_block" "NodeName")
    local state=$(extract_value "$node_block" "State")
    local cpu_total=$(extract_value "$node_block" "CPUTot")
    local cpu_alloc=$(extract_value "$node_block" "CPUAlloc")
    local mem_total=$(extract_value "$node_block" "RealMemory")
    local mem_alloc=$(extract_value "$node_block" "AllocMem")
    local partitions=$(extract_value "$node_block" "Partitions")
    local features=$(extract_value "$node_block" "ActiveFeatures")
    local gres=$(extract_value "$node_block" "Gres")
    local alloc_tres=$(extract_value "$node_block" "AllocTRES")
    
    # Skip if not a GPU node
    if [[ "$gres" != *"gpu"* ]]; then
        return
    fi
    
    # Handle memory units (remove G suffix if present)
    mem_alloc=${mem_alloc%G}
    
    # Calculate derived metrics
    local cpu_free=$((cpu_total - cpu_alloc))
    local cpu_util=$(awk "BEGIN {printf \"%.1f\", $cpu_alloc * 100 / $cpu_total}")
    
    local mem_free=$((mem_total - mem_alloc))
    local mem_util=$(awk "BEGIN {printf \"%.1f\", $mem_alloc * 100 / $mem_total}")
    
    # Safe GPU parsing section in process_single_node function
    # Extract GPU information with error handling
    local gpu_info=$(extract_gpu_info "$gres")
    if [[ -z "$gpu_info" ]]; then
        gpu_info="none,0"
    fi

    IFS=',' read -r gpu_type gpu_total <<< "$gpu_info"

    # Ensure variables have valid values
    if [[ -z "$gpu_type" ]]; then
        gpu_type="none"
    fi
    if [[ -z "$gpu_total" || ! "$gpu_total" =~ ^[0-9]+$ ]]; then
        gpu_total=0
    fi

    local gpu_alloc=$(extract_alloc_gpus "$alloc_tres")
    if [[ -z "$gpu_alloc" || ! "$gpu_alloc" =~ ^[0-9]+$ ]]; then
        gpu_alloc=0
    fi

    local gpu_free=$((gpu_total - gpu_alloc))
    local gpu_util=0
    if [[ $gpu_total -gt 0 ]]; then
        gpu_util=$(awk "BEGIN {printf \"%.1f\", $gpu_alloc * 100 / $gpu_total}")
    fi
    
    # Output in requested format
    if [[ "$output_format" == "csv" ]]; then
        echo "$node_name,$state,$cpu_total,$cpu_alloc,$cpu_free,$cpu_util,$mem_total,$mem_alloc,$mem_free,$mem_util,$gpu_type,$gpu_total,$gpu_alloc,$gpu_free,$gpu_util,$partitions,$features"
    elif [[ "$output_format" == "json" ]]; then
        if [[ "$first_node" == "false" ]]; then
            echo ","
        fi
        cat << EOF
  {
    "node_name": "$node_name",
    "state": "$state",
    "cpu": {
      "total": $cpu_total,
      "allocated": $cpu_alloc,
      "free": $cpu_free,
      "utilization": $cpu_util
    },
    "memory": {
      "total": $mem_total,
      "allocated": $mem_alloc,
      "free": $mem_free,
      "utilization": $mem_util
    },
    "gpu": {
      "type": "$gpu_type",
      "total": $gpu_total,
      "allocated": $gpu_alloc,
      "free": $gpu_free,
      "utilization": $gpu_util
    },
    "partitions": "$partitions",
    "features": "$features"
  }
EOF
    else
        # Default tabular format
        printf "%-15s %-10s %2d/%-2d (%5.1f%%) %5d/%-6d (%5.1f%%) %-6s %1d/%-1d (%5.1f%%) %s\n" \
            "$node_name" "$state" "$cpu_alloc" "$cpu_total" "$cpu_util" \
            "$mem_alloc" "$mem_total" "$mem_util" "$gpu_type" \
            "$gpu_alloc" "$gpu_total" "$gpu_util" "$partitions"
    fi
}

show_summary() {
    echo "=== GPU Node Resource Summary ==="
    echo
    
    local total_nodes=0
    local idle_nodes=0
    local mixed_nodes=0
    local alloc_nodes=0
    local down_nodes=0
    local other_nodes=0
    local total_gpus=0
    local alloc_gpus=0
    
    # Process each node for summary
    local current_node=""
    local in_node=false
    
    while IFS= read -r line; do
        if [[ "$line" =~ ^NodeName= ]]; then
            if [[ -n "$current_node" ]]; then
                # FIXED: Call with correct number of parameters
                process_node_summary_fixed "$current_node" total_nodes idle_nodes mixed_nodes alloc_nodes down_nodes other_nodes total_gpus alloc_gpus
            fi
            current_node="$line"
            in_node=true
        elif [[ "$in_node" == true ]]; then
            current_node="$current_node $line"
        fi
    done < <(get_scontrol_output)
    
    # Process the last node
    if [[ -n "$current_node" ]]; then
        process_node_summary_fixed "$current_node" total_nodes idle_nodes mixed_nodes alloc_nodes down_nodes other_nodes total_gpus alloc_gpus
    fi
    
    echo "GPU Nodes by State:"
    echo "  Total: $total_nodes"
    echo "  Idle: $idle_nodes"
    echo "  Mixed: $mixed_nodes"
    echo "  Allocated: $alloc_nodes"
    echo "  Down: $down_nodes"
    echo "  Other: $other_nodes"
    echo
    
    local free_gpus=$((total_gpus - alloc_gpus))
    local gpu_utilization=0
    if [[ $total_gpus -gt 0 ]]; then
        gpu_utilization=$(awk "BEGIN {printf \"%.1f\", $alloc_gpus * 100 / $total_gpus}")
    fi
    
    echo "GPU Resources:"
    echo "  Total GPUs: $total_gpus"
    echo "  Allocated GPUs: $alloc_gpus"
    echo "  Free GPUs: $free_gpus"
    echo "  Utilization: ${gpu_utilization}%"
    echo
}

# Add this NEW function to your script (don't replace the existing process_node_summary)
process_node_summary_fixed() {
    local node_block="$1"
    local -n total_ref=$2
    local -n idle_ref=$3
    local -n mixed_ref=$4
    local -n alloc_ref=$5
    local -n down_ref=$6
    local -n other_ref=$7
    local -n total_gpu_ref=$8
    local -n alloc_gpu_ref=$9
    
    local gres=$(extract_value "$node_block" "Gres")
    [[ "$gres" != *"gpu"* ]] && return
    
    local state=$(extract_value "$node_block" "State")
    local alloc_tres=$(extract_value "$node_block" "AllocTRES")
    
    ((total_ref++))
    
    # Use improved state categorization
    local state_category=$(get_node_state_category "$state")
    case "$state_category" in
        "IDLE") ((idle_ref++)) ;;
        "MIXED") ((mixed_ref++)) ;;
        "ALLOCATED") ((alloc_ref++)) ;;
        "DOWN") ((down_ref++)) ;;
        *) ((other_ref++)) ;;
    esac
    
    # FIXED: Safe GPU info extraction without IFS read
    local gpu_info=$(extract_gpu_info "$gres")
    
    # Ensure we have a valid gpu_info string
    if [[ -z "$gpu_info" || "$gpu_info" == "," ]]; then
        gpu_info="none,0"
    fi
    
    # FIXED: Safe parsing using parameter expansion
    local gpu_type="none"
    local gpu_total=0
    
    # Parse using parameter expansion instead of IFS read
    if [[ "$gpu_info" == *","* ]]; then
        gpu_type="${gpu_info%,*}"  # Everything before the last comma
        gpu_total="${gpu_info##*,}" # Everything after the last comma
    fi
    
    # Validate extracted values
    if [[ -z "$gpu_type" ]]; then
        gpu_type="none"
    fi
    if [[ ! "$gpu_total" =~ ^[0-9]+$ ]]; then
        gpu_total=0
    fi
    
    local gpu_alloc=$(extract_alloc_gpus "$alloc_tres")
    if [[ ! "$gpu_alloc" =~ ^[0-9]+$ ]]; then
        gpu_alloc=0
    fi
    
    ((total_gpu_ref += gpu_total))
    ((alloc_gpu_ref += gpu_alloc))
}

# Improved state handling
get_node_state_category() {
    local state="$1"
    
    # Handle compound states by checking primary state
    case "$state" in
        IDLE*) echo "IDLE" ;;
        MIXED*) echo "MIXED" ;;
        ALLOCATED*) echo "ALLOCATED" ;;
        DOWN*|DRAIN*|*DRAIN*|*NOT_RESPONDING*) echo "DOWN" ;;
        COMPLETING*) echo "COMPLETING" ;;
        *) echo "OTHER" ;;
    esac
}

# Replace the existing process_node_summary function with this fixed version
process_node_summary() {
    local node_block="$1"
    local -n total_ref=$2
    local -n idle_ref=$3
    local -n mixed_ref=$4
    local -n alloc_ref=$5
    local -n down_ref=$6
    local -n other_ref=$7
    local -n total_gpu_ref=$8
    local -n alloc_gpu_ref=$9
    
    local gres=$(extract_value "$node_block" "Gres")
    [[ "$gres" != *"gpu"* ]] && return
    
    local state=$(extract_value "$node_block" "State")
    local alloc_tres=$(extract_value "$node_block" "AllocTRES")
    
    ((total_ref++))
    
    # Use improved state categorization
    local state_category=$(get_node_state_category "$state")
    case "$state_category" in
        "IDLE") ((idle_ref++)) ;;
        "MIXED") ((mixed_ref++)) ;;
        "ALLOCATED") ((alloc_ref++)) ;;
        "DOWN") ((down_ref++)) ;;
        *) ((other_ref++)) ;;
    esac
    
    # FIXED: Safe GPU info extraction
    local gpu_info=$(extract_gpu_info "$gres")
    
    # Ensure we have a valid gpu_info string
    if [[ -z "$gpu_info" || "$gpu_info" == "," ]]; then
        gpu_info="none,0"
    fi
    
    # FIXED: Safe parsing instead of IFS read
    local gpu_type=""
    local gpu_total=0
    
    # Parse using parameter expansion instead of IFS read
    if [[ "$gpu_info" == *","* ]]; then
        gpu_type="${gpu_info%,*}"  # Everything before the last comma
        gpu_total="${gpu_info##*,}" # Everything after the last comma
    else
        gpu_type="none"
        gpu_total=0
    fi
    
    # Validate extracted values
    if [[ -z "$gpu_type" ]]; then
        gpu_type="none"
    fi
    if [[ ! "$gpu_total" =~ ^[0-9]+$ ]]; then
        gpu_total=0
    fi
    
    local gpu_alloc=$(extract_alloc_gpus "$alloc_tres")
    if [[ ! "$gpu_alloc" =~ ^[0-9]+$ ]]; then
        gpu_alloc=0
    fi
    
    ((total_gpu_ref += gpu_total))
    ((alloc_gpu_ref += gpu_alloc))
}

# Function to recommend optimal job size
recommend_job_size() {
    echo "=== Job Size Recommendations ==="
    echo
    
    local idle_gpus=0
    local mixed_free_gpus=0
    local current_node=""
    local in_node=false
    
    while IFS= read -r line; do
        if [[ "$line" =~ ^NodeName= ]]; then
            if [[ -n "$current_node" ]]; then
                process_node_recommendation "$current_node" idle_gpus mixed_free_gpus
            fi
            current_node="$line"
            in_node=true
        elif [[ "$in_node" == true ]]; then
            current_node="$current_node $line"
        fi
    done < <(get_scontrol_output)
    
    # Process the last node
    if [[ -n "$current_node" ]]; then
        process_node_recommendation "$current_node" idle_gpus mixed_free_gpus
    fi
    
    local total_available=$((idle_gpus + mixed_free_gpus))
    
    echo "Available Resources:"
    echo "  GPUs on idle nodes: $idle_gpus"
    echo "  GPUs free on mixed nodes: $mixed_free_gpus"
    echo "  Total available GPUs: $total_available"
    echo
    echo "Recommendations:"
    if [[ $total_available -gt 0 ]]; then
        echo "  Small job (1-2 GPUs): Immediate availability"
        if [[ $total_available -ge 4 ]]; then
            echo "  Medium job (4 GPUs): Good availability"
        fi
        if [[ $total_available -ge 8 ]]; then
            echo "  Large job (8+ GPUs): Consider node distribution"
        fi
    else
        echo "  No GPUs currently available. Consider scheduling for later."
    fi
}

# Add this debug function to see all node states:
debug_states() {
    echo "=== All GPU Node States ==="
    echo
    
    local current_node=""
    local in_node=false
    
    while IFS= read -r line; do
        if [[ "$line" =~ ^NodeName= ]]; then
            if [[ -n "$current_node" ]]; then
                local gres=$(extract_value "$current_node" "Gres")
                if [[ "$gres" == *"gpu"* ]]; then
                    local node_name=$(extract_value "$current_node" "NodeName")
                    local state=$(extract_value "$current_node" "State")
                    printf "%-15s %s\n" "$node_name" "$state"
                fi
            fi
            current_node="$line"
            in_node=true
        elif [[ "$in_node" == true ]]; then
            current_node="$current_node $line"
        fi
    done < <(get_scontrol_output)
    
    # Process the last node
    if [[ -n "$current_node" ]]; then
        local gres=$(extract_value "$current_node" "Gres")
        if [[ "$gres" == *"gpu"* ]]; then
            local node_name=$(extract_value "$current_node" "NodeName")
            local state=$(extract_value "$current_node" "State")
            printf "%-15s %s\n" "$node_name" "$state"
        fi
    fi
}

# Helper function for recommendation processing
process_node_recommendation() {
    local node_block="$1"
    local -n idle_ref=$2
    local -n mixed_ref=$3
    
    local gres=$(extract_value "$node_block" "Gres")
    [[ "$gres" != *"gpu"* ]] && return
    
    local state=$(extract_value "$node_block" "State")
    local alloc_tres=$(extract_value "$node_block" "AllocTRES")
    
    local gpu_info=$(extract_gpu_info "$gres")
    IFS=',' read -r gpu_type gpu_total <<< "$gpu_info"
    local gpu_alloc=$(extract_alloc_gpus "$alloc_tres")
    local gpu_free=$((gpu_total - gpu_alloc))
    
    case "$state" in
        "IDLE") ((idle_ref += gpu_total)) ;;
        "MIXED") ((mixed_ref += gpu_free)) ;;
    esac
}

# Main script logic
case "${1:-table}" in
    "csv")
        process_nodes "csv"
        ;;
    "json")
        process_nodes "json"
        ;;
    "summary")
        show_summary
        ;;
    "recommend")
        recommend_job_size
        ;;
    "table"|*)
        echo "=== GPU Node Status ==="
        printf "%-15s %-10s %-15s %-18s %-8s %-13s %s\n" \
            "Node" "State" "CPU (Used/Total)" "Memory (Used/Total)" "GPU Type" "GPU (Used/Total)" "Partitions"
        echo "$(printf '%.0s-' {1..100})"
        process_nodes "table"
        echo
        echo "Usage: $0 [table|csv|json|summary|recommend]"
        ;;
    "debug")
        debug_states
        ;;
esac





