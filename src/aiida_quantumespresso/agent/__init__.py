"""aiida-agents plugin: what aiida-quantumespresso contributes to the agent.

aiida-agents discovers this through the ``aiida_agents.plugins`` entry point
(see ``pyproject.toml``), which resolves to this module. Its :data:`name` and
the :func:`rag_corpora` / :func:`tools` / :func:`prompt_fragment` functions are
the provider; the module *is* the plugin, since a provider holds no state and
discovery only reads these attributes. Every hook is optional and read
defensively.

This is the only place in aiida-quantumespresso allowed to import
``aiida_agents``: it loads only during aiida-agents' discovery, so aiida-agents
stays an optional extra, never a runtime dependency of the core.

Real: :func:`rag_corpora`, :func:`prompt_fragment`. Stubbed with a TODO:
:func:`tools`.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from aiida_agents.plugins import AgentTool, RagCorpus

__all__ = ['name', 'prompt_fragment', 'rag_corpora', 'tools']

#: Name aiida-agents shows for this plugin's tools and doc citations.
name = 'quantumespresso'

#: Distribution whose version stamps the RAG corpus identity and its docs ref.
_DISTRIBUTION = 'aiida-quantumespresso'


def _distribution_version() -> str:
    """Return the installed ``aiida-quantumespresso`` version.

    The corpus identity (see :class:`~aiida_agents.plugins.RagCorpus`) is keyed
    on this, so a version bump resolves to a different collection and triggers a
    rebuild rather than silently disagreeing with the installed code.

    :return: The installed distribution version, or the in-tree ``__version__``
        when the distribution metadata is unavailable (a bare source checkout).
    """
    try:
        return version(_DISTRIBUTION)
    except PackageNotFoundError:
        from aiida_quantumespresso import __version__

        return __version__


def rag_corpora() -> list[RagCorpus]:
    """Return the documentation corpora to index and cite.

    Registers the published Quantum ESPRESSO plugin documentation, versioned
    from the installed distribution, cloned at the matching release tag, and
    linked back to Read the Docs so retrieved passages cite an openable URL.

    :return: A single-element list with the Quantum ESPRESSO docs corpus.
    """
    release = _distribution_version()
    return [
        RagCorpus(
            name='quantumespresso',
            version=release,
            docs_repo='https://github.com/aiidateam/aiida-quantumespresso.git',
            # The tag the corpus is rendered from. A tagged release resolves
            # cleanly; a development install (e.g. '5.0.1.dev0') has no matching
            # tag, so pin an explicit ref here when iterating from an untagged
            # checkout.
            docs_ref=f'v{release}',
            docs_subdir='docs',
            # Where the same docs are published, so a hit is cited with a link.
            # '{version}' is filled from docs_ref and '{page}' from the page
            # path, so the linked page is the one the corpus was rendered from.
            docs_url='https://aiida-quantumespresso.readthedocs.io/en/{version}/{page}.html',
        )
    ]


def tools() -> list[AgentTool]:
    """Return the tools this plugin contributes to the agent.

    :return: Currently empty.

    .. todo::

        Contribute Quantum ESPRESSO specific read tools that parse QE's own
        output, which the code-agnostic core tool layer cannot. The obvious
        first one reads a pw.x SCF convergence trace (accuracy series, trend,
        stall) so the agent can tell a run that merely ran out of iterations
        from one that oscillated. A tool that changes state (submitting,
        storing) must be declared with ``AgentTool(fn=..., writes=True)`` so it
        is registered behind the human-in-the-loop approval gate.
    """
    return []


#: Quantum ESPRESSO guidance appended to the agent's system prompt under a
#: '### quantumespresso' heading. Conventions and physics the code-agnostic core
#: cannot know; well within ``discovery.MAX_PROMPT_FRAGMENT_CHARS`` (2000).
_PROMPT_FRAGMENT = """\
Quantum ESPRESSO here means the pw.x code, driven through aiida-quantumespresso \
work chains rather than bare calculations, so restarts and error handling apply.

Choosing a work chain:
- PwBaseWorkChain: a single pw.x run (scf, relax, vc-relax, bands or nscf) with \
automatic restart on the common failures.
- PwRelaxWorkChain: geometry optimisation, re-running vc-relax until the cell is \
self-consistent.
- PwBandsWorkChain: scf then a bands run (optionally preceded by a relax).
- PdosWorkChain: scf, nscf, then projwfc for a projected density of states.

Setting up inputs:
- Prefer get_builder_from_protocol(code, structure, protocol=...), with a \
protocol of 'fast', 'moderate' or 'precise', over setting parameters by hand; \
it also selects pseudopotentials and cutoffs.
- Pseudopotentials come from an aiida-pseudo family (SSSP or PseudoDojo). Take \
ecutwfc and ecutrho from that family's recommended cutoffs for the elements \
present; do not invent them.

When SCF will not converge it is usually mixing or smearing, not too few steps: \
lower mixing_beta (around 0.2 to 0.3), and for metals use smearing occupations \
with a small degauss (about 0.01 to 0.02 Ry). Raising electron_maxstep alone \
rarely fixes an oscillating run.
"""


def prompt_fragment() -> str:
    """Return Quantum ESPRESSO guidance appended to the agent's system prompt.

    Conventions and physics only this plugin knows: which work chain to run,
    that protocols and pseudo families supply cutoffs, and the usual SCF
    convergence remedies. Kept well within ``discovery.MAX_PROMPT_FRAGMENT_CHARS``.

    :return: The guidance fragment.
    """
    return _PROMPT_FRAGMENT
