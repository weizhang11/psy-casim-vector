
from psyclone.psyGen import Transformation, CodedKern
from psyclone.psyir.nodes import Routine, Loop
from psyclone.transformations import ACCRoutineTrans, ACCLoopTrans, KernelImportsToArguments
from utils import enhance_tree_information, normalise_loops
from utils_casim import ACCRoutineTrans_Force

def trans(psy):
    """Apply `acc routine` to allow OpenACC loops."""
    
    acc_routine_trans = ACCRoutineTrans()
    acc_loop_trans = ACCLoopTrans()
    acc_routine_trans_force = ACCRoutineTrans_Force()


    for routine in psy.walk(Routine):  # Get all subroutines

        enhance_tree_information(routine)

        normalise_loops(
                routine,
                hoist_local_arrays=True,
                convert_array_notation=True,
                convert_range_loops=True,
                hoist_expressions=True
        )

        #acc_routine_trans.apply(routine, {"force": True, "parallelism": "vector"})  # Apply `!$acc routine`
        #print(f"{routine.name} as GPU-enabled with 'acc routine vector'")

        acc_routine_trans_force.apply(routine, {"force": True, "parallelism": "vector"})  # Apply `!$acc routine`
        print(f"Forced {routine.name} as GPU-enabled with 'acc routine vector'")
        
        for loop in routine.walk(Loop):  # Apply loop transformation
            acc_loop_trans.apply(loop, {"vector": True, "independent": False})  

    return psy

