Workchain development
+++++++++++++++++++++

The ``BaseRestartWorkChain``
============================

During the development of the workchains for the ``pw.x`` and ``ph.x`` codes, with their respective |PhCalculation| and
|PwCalculation|, a common pattern emerged. Calculations for both codes have the concept of a restart. If convergence at
the end of a calculation was not reached, a next calculation can restart from the previous calculation and try again to
achieve convergence. Of course calculations could also fail for a multitude of other reasons that would need to be dealt
with. For example the calculation could fail to submit to the scheduler, or the calculation could fail for a variety of
reasons. For each launched calculation, the workchain needed to check for all of these cases and act appropriately. To
prevent from having to implement this same logic, it was abstracted to the |BaseRestartWorkChain|.

Even though the |BaseRestartWorkChain| subclasses the :py:class:`~aiida.engine.processes.workchains.WorkChain` class, it is
technically in itself not a runnable workchain, for example because it does not define any inputs or an outline in its
``spec``. Rather, this class should be subclassed by actual workchains, such as the |PhBaseWorkChain| and the
|PwBaseWorkChain|, which can then leverage its predefined outline methods.

There are just a few mandatory steps necessary to make this work.

Define the ``_calculation_class``
---------------------------------

A workchain implementation that subclasses the |BaseRestartWorkChain| needs to override the ``_process_class`` attribute with the correct calculation class, e.g. |PhCalculation| or
|PwCalculation|. The workchain will submit instances of this job calculation class.

Define the ``spec.outline``
---------------------------
The outline of the workchain spec needs to be defined. The most basic and minimally required example would look like
the following snippet::

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.outline(
            cls.setup,
            while_(cls.should_run_calculation)(
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )

All these five outline class methods have already been implemented by the |BaseRestartWorkChain| and one just needs to
add them to the outline. Note that ``cls.setup`` **has** to be called and if it is overridden, one needs to make sure
to make the super call or the workchain will break.

The developer may then of course add methods to the outline or override the class methods of
the |BaseRestartWorkChain|. Take a look at for example the |PwBaseWorkChain| to see how additional outline methods are
added, defined and used.

Define the ``self.ctx.inputs``
------------------------------
For the ``run_calculation`` method to work, the user has to define a dictionary of inputs that are supposed to be
passed to the calculation that it will submit. This dictionary has to be defined in the context under the key ``inputs``.
For example in an outline method ``prepare_calculation`` one could do the following::

    def prepare_calculation(self):
        """
        Prepare the inputs dictionary for the run_calculation call
        """
        self.ctx.inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.parameters
        }

The values in the dictionary need to be instances of :class:`aiida.orm.nodes.data.Data` or a dictionary. The |BaseRestartWorkChain|
will make sure that bare dictionaries will be wrapped in :class:`~aiida.orm.nodes.data.dict.Dict` instances.

Error handling
==============
In the previous paragraph we explained how the |BaseRestartWorkChain| can be used to quickly write a workchain that will
automatically deal with the most basic functionality of submitting a calculation and dealing with any generic failures
such as submission failures. All of this is agnostic of the actual calculation, but of course the handling of actual
calculation failures is going to be calculation specific. To allow the user to register error handling functions to the
workchain, the |register_error_handler| decorator is defined, which can promote a plain python function to an error
handling class method for any workchain that extends the |BaseRestartWorkChain|.

Defining an error handler
-------------------------
To define an error handler for a workchain, one should write a function that takes two arguments: ``self`` and
``calculation``. The first will refer to the instance of the workchain and the second to the instance of the calculation
that failed. The body should typically be a single conditional that matches some specific potential calculation failure.
If it is matched, the handler can change the inputs to fix the problem and or report some messages. A minimal example
would look something like this::

    def _handler_error_generic(self, calculation):
        if 'typical error encountered' in calculation.res.warnings:
            self.ctx.inputs.parameters['alpha'] = 2.0
            self.report('incorrect value for alpha, reset it to 2.0')
            return ErrorHandlerReport(True, True)

If the conditional is matched, the inputs dictionary in the context is updated and we fire a report so it is logged.
Finally a ``ProcessHandlerReport`` is returned to tell the |BaseRestartWorkChain| that the error was handled and no
further error handlers should be called and the next iteration should be performed. If the ``calculation`` can be
restarted from in the next iteration, despite the calculation failure, one can set it to the ``restart_calc`` member of
the context. This will cause the workchain to automatically use this calculation to restart from::

    def _handler_error_generic(self, calculation):
        self.ctx.restart_calc = calcuation

Now how do we add this error handler to the actual workchain?

The ``register_error_handler`` decorator
----------------------------------------
To add an error handling function to a particular workchain class, one should use the |register_error_handler|
decorator. In the same file were the workchain in question is defined, one can write something like the following::

    @register_error_handler(PhBaseWorkChain, 300)
    def _handle_error_exceeded_maximum_walltime(self, calculation):
        """
        Calculation ended nominally but ran out of allotted wall time
        """
        if 'Maximum CPU time exceeded' in calculation.res.warnings:
            self.ctx.restart_calc = calculation
            self.report('PhCalculation<{}> exceeded max wall time, restarting'
                .format(calculation.pk))
            return ErrorHandlerReport(True, True)

The decorator takes two arguments: the workchain class to which the handler should be added and an integer indicating
the priority with which it should be called with respect to other handlers. This allows the user to control the order
in which handlers will be called. Handlers with a higher priority will be called first.
That is all. The decorator will make sure that the workchain class gets the function as a class method and in the
``inspect_calculation`` call, when a calculation has failed, the workchain will loop over
all the registered error handlers and call them.

The ``_error_handler_entry_point``
----------------------------------
In the previous paragraph, we explained how the |register_error_handler| decorator could register a function as an
error handler for a |BaseRestartWorkChain|. One condition was that the function was defined in the same file as the
workchain class itself. This is because the decorator, and therefore the registration, only gets performed when the
function is imported. Putting it in the same file as the workchain class guarantees that this happens. But what if we
do not have write access to that file?

To solve this problem, the |BaseRestartWorkChain| has the |error_handler_entry_point| attribute. The subclassing workchain
can define an entry point category, for example::

    _error_handler_entry_point = 'aiida_quantumespresso.workflow_error_handlers.pw.base'

One can then register entry points to this category that point to a file, in which additional error handler are defined
with the |register_error_handler| handler. Upon construction of the workchain, the ``aiida.common.pluginloader`` will be
used to import the files registered under that entry point, causing the decorators to be called and the error handlers
to be registered with the workchain.

To add entries to the error handler category from another package, simply define it in the ``setup.json``::

    "entry_points": {
        "aiida_quantumespresso.workflow_error_handlers.pw.base": [
            "epfl = aiida_quantumespresso_epfl.workflows.pw.base"
        ]
    }

where the ``aiida_quantumespresso_epfl.workflows.pw.base`` file contains the additional decorated error handlers.

.. |error_handler_entry_point| replace:: `_error_handler_entry_point`
.. |register_error_handler| replace:: :py:func:`~aiida.engine.processes.workchains.utils.process_handler`
.. |BaseRestartWorkChain| replace:: :py:class:`~aiida.engine.processes.workchains.restart.BaseRestartWorkChain`
.. |PhCalculation| replace:: :py:class:`~aiida_quantumespresso.calculations.ph.PhCalculation`
.. |PwCalculation| replace:: :py:class:`~aiida_quantumespresso.calculations.pw.PwCalculation`
.. |PhBaseWorkChain| replace:: :py:class:`~aiida_quantumespresso.workflows.ph.base.PhBaseWorkChain`
.. |PwBaseWorkChain| replace:: :py:class:`~aiida_quantumespresso.workflows.pw.base.PwBaseWorkChain`
