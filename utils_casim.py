import abc

from psyclone import psyGen
from psyclone.configuration import Config
from psyclone.core import Signature, VariablesAccessInfo
from psyclone.domain.lfric import (KernCallArgList, LFRicConstants,
                                   LFRicInvokeSchedule, LFRicKern, LFRicLoop)
from psyclone.dynamo0p3 import LFRicHaloExchangeEnd, LFRicHaloExchangeStart
from psyclone.errors import InternalError
from psyclone.gocean1p0 import GOInvokeSchedule
from psyclone.psyGen import (Transformation, CodedKern, Kern, InvokeSchedule,
                             BuiltIn)
from psyclone.psyir.nodes import (
    ACCDataDirective, ACCDirective, ACCEnterDataDirective, ACCKernelsDirective,
    ACCLoopDirective, ACCParallelDirective, ACCRoutineDirective,
    Call, CodeBlock, Directive, Loop, Node,
    OMPDeclareTargetDirective, OMPDirective, OMPMasterDirective,
    OMPParallelDirective, OMPParallelDoDirective, OMPSerialDirective,
    OMPSingleDirective, OMPTaskloopDirective, PSyDataNode, Reference,
    Return, Routine, Schedule)
from psyclone.psyir.nodes.array_mixin import ArrayMixin
from psyclone.psyir.nodes.structure_member import StructureMember
from psyclone.psyir.nodes.structure_reference import StructureReference
from psyclone.psyir.symbols import (
    ArgumentInterface, DataSymbol, INTEGER_TYPE, ScalarType, Symbol,
    SymbolError, UnresolvedType)
from psyclone.psyir.transformations.loop_trans import LoopTrans
from psyclone.psyir.transformations.omp_loop_trans import OMPLoopTrans
from psyclone.psyir.transformations.parallel_loop_trans import (
    ParallelLoopTrans)
from psyclone.psyir.transformations.region_trans import RegionTrans
from psyclone.psyir.transformations.transformation_error import (
    TransformationError)

class ACCRoutineTrans_Force(Transformation):
    '''
    Transform a kernel or routine by adding a "!$acc routine" directive
    (causing it to be compiled for the OpenACC accelerator device).
    For example:

    >>> from psyclone.parse.algorithm import parse
    >>> from psyclone.psyGen import PSyFactory
    >>> api = "gocean"
    >>> ast, invokeInfo = parse(GOCEAN_SOURCE_FILE, api=api)
    >>> psy = PSyFactory(api).create(invokeInfo)
    >>>
    >>> from psyclone.transformations import ACCRoutineTrans
    >>> rtrans = ACCRoutineTrans()
    >>>
    >>> schedule = psy.invokes.get('invoke_0').schedule
    >>> # Uncomment the following line to see a text view of the schedule
    >>> # print(schedule.view())
    >>> kern = schedule.children[0].children[0].children[0]
    >>> # Transform the kernel
    >>> rtrans.apply(kern)

    '''
    def apply(self, node, options=None):
        '''
        Add the '!$acc routine' OpenACC directive into the code of the
        supplied Kernel (in a PSyKAl API such as GOcean or LFRic) or directly
        in the supplied Routine.

        :param node: the kernel call or routine implementation to transform.
        :type node: :py:class:`psyclone.psyGen.Kern` |
                    :py:class:`psyclone.psyir.nodes.Routine`
        :param options: a dictionary with options for transformations.
        :type options: Optional[Dict[str, Any]]
        :param bool options["force"]: whether to allow routines with
            CodeBlocks to run on the GPU.
        :param str options["parallelism"]: the level of parallelism that the
            target routine (or a callee) exposes. One of "seq" (the default),
            "vector", "worker" or "gang".

        '''
        # Check that we can safely apply this transformation
        self.validate(node, options)

        if isinstance(node, Kern):
            # Flag that the kernel has been modified
            node.modified = True

            # Get the schedule representing the kernel subroutine
            routine = node.get_kernel_schedule()
        else:
            routine = node

        # Insert the directive to the routine if it doesn't already exist
        for child in routine.children:
            if isinstance(child, ACCRoutineDirective):
                return  # The routine is already marked with ACCRoutine

        para = options.get("parallelism", "seq") if options else "seq"

        routine.children.insert(0, ACCRoutineDirective(parallelism=para))

#    def validate(self, node, options=None):
#        '''
#        Perform checks that the supplied kernel or routine can be transformed.
#
#        :param node: the kernel or routine which is the target of this
#            transformation.
#        :type node: :py:class:`psyclone.psyGen.Kern` |
#                    :py:class:`psyclone.psyir.nodes.Routine`
#        :param options: a dictionary with options for transformations.
#        :type options: Optional[Dict[str, Any]]
#        :param bool options["force"]: whether to allow routines with
#            CodeBlocks to run on the GPU.
#
#        :raises TransformationError: if the node is not a kernel or a routine.
#        :raises TransformationError: if the target is a built-in kernel.
#        :raises TransformationError: if it is a kernel but without an
#            associated PSyIR.
#        :raises TransformationError: if any of the symbols in the kernel are
#            accessed via a module use statement.
#        :raises TransformationError: if the kernel contains any calls to other
#            routines.
#        :raises TransformationError: if the 'parallelism' option is supplied
#            but is not a recognised level of parallelism.
#
#        '''
#
#        force = options.get("force", False) if options else False
#
#        if not force:
#           super().validate(node, options)
#           self.validate_it_can_run_on_gpu(node, options)
#
#        if options and "parallelism" in options:
#            para = options["parallelism"]
#            if para not in ACCRoutineDirective.SUPPORTED_PARALLELISM:
#                raise TransformationError(
#                    f"{self.name}: '{para}' is not a supported level of "
#                    f"parallelism. Should be one of "
#                    f"{ACCRoutineDirective.SUPPORTED_PARALLELISM}")

__all__ = [
   "ACCRoutineTrans_Force",
]
