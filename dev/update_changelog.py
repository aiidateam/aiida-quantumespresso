#!/bin/bash
# -*- coding: utf-8 -*-
"""Script for automatically updating the `CHANGELOG.md` based on the commits since the latest release tag."""
from pathlib import Path
import re
import subprocess

DEFAULT_CHANGELOG_SECTIONS = """
### ‼️ Breaking changes


### ✨ New features


### 🗑️ Deprecation


### 👌 Improvements


### 🐛 Bug Fixes


### 📚 Documentation


### 🔧 Maintenance


### ⬆️ Update dependencies


### ♻️ Refactor

"""


def update_changelog():
    """Update the `CHANGELOG.md` for a first draft of the release."""

    print('🔍 Checking the current version number')
    with Path('CHANGELOG.md').open('r', encoding='utf8') as handle:
        current_changelog = handle.read()

    with Path('src/aiida_quantumespresso/__init__.py').open('r', encoding='utf8') as handle:
        version = re.search(r"__version__ = '(?P<version_number>\d+.+)'", handle.read()).groupdict()['version_number']

    if str(version) in current_changelog:
        print('🛑 Current version already in `CHANGELOG.md`. Skipping...')
        return

    print('⬆️ Found updated version number, adapting `CHANGELOG.md`.')
    tags = subprocess.run(['git', 'tag'], capture_output=True, check=False).stdout

    tag_pattern = re.compile(r'(v\d\.\d\.\d)\n')
    tags = tags.decode()

    latest_tag = tag_pattern.findall(tags)[-1]

    print(f'🔄 Comparing with latest tag `{latest_tag}`.')
    commits = subprocess.run(['git', 'log', "--pretty=format:'%s %h %H'", f'{latest_tag}..origin/main'],
                             capture_output=True,
                             check=False).stdout
    commits = commits.decode()

    pr_pattern = re.compile(r'\(\S(?P<pr_number>\d+)\)')

    changelog_commits = []

    for commit in commits.splitlines():

        # Remove the PR number from the commit message
        pr_match = pr_pattern.search(commit)

        if pr_match is not None:
            pr_number = pr_match.groupdict()['pr_number']
            commit = commit.replace(fr'(#{pr_number})', '')

        # Add the commit hash (short) to link to the
        # 0aba276f7042d51b91d5699530de7336cd62e2c2
        commit = commit.split()
        hash_long = commit.pop()
        hash_short = commit.pop()
        commit.append(f'[[{hash_short}](https://github.com/aiidateam/aiida-quantumespresso/commit/{hash_long})]')
        commit = commit.strip("'")
        changelog_commits.append(' '.join(commit))

    changelog_message = f'## v{version}\n' + DEFAULT_CHANGELOG_SECTIONS

    for commit_title in changelog_commits:
        changelog_message += f'\n* {commit_title}'

    with Path('CHANGELOG.md').open('w', encoding='utf8') as handle:
        handle.write(changelog_message + '\n\n' + current_changelog)

    print("🚀 Success! Finalise the `CHANGELOG.md` and let's get this baby released.")


if __name__ == '__main__':
    update_changelog()
