# PolyzyMD Docs Collaboration Guide

This directory contains operating instructions for documentation collaborators.

Current collaborators:

- `docs/agents/diataxis-doc-writer.md`
- `docs/agents/docs-ia-curator.md`

## Mission

Overhaul the PolyzyMD documentation so users and LLMs can quickly find the
right information, with the docs organized around Diataxis and shaped by good
scientific-software documentation practice.

Primary references:

- Diataxis: `https://diataxis.fr/`
- OMSF Documentation Playbook: `https://playbooks.omsf.io/documentation/`
- Style references:
  - `https://docs.openfree.energy/en/latest/tutorials/ahfe_tutorial.html`
  - `https://docs.openforcefield.org/en/latest/install.html`

## Target Information Architecture

PolyzyMD should move toward this top-level structure:

- `Home`
- `Get Started`
- `Tutorials`
- `How-To Guides`
- `Reference`
- `Explanation`
- `Contributor Guide`

Contributor documentation stays top-level, not hidden.

## Project-Specific Content Rules

- Keep `pixi` as the primary package manager and installation story.
- Keep stable vs experimental analysis workflows clearly distinguished.
- Experimental metrics remain documented, but must be explicitly labeled as
  experimental and not presented as settled default science.
- Do not invent workflows, commands, or capabilities not present in the repo.
- Prefer concise, high-signal landing pages over long link dumps.
- Keep API reference pages in `Reference`, not mixed into tutorials.

## Current Workflow

This overhaul happens in two phases:

1. `Planning and classification`
   - collaborators inspect the current docs
   - collaborators propose IA, page moves, splits, and rewrites
   - the main agent summarizes the plan for the user
   - no major docs rewrite begins until the user approves

2. `Implementation`
   - once approved, collaborators help execute the rewrite in batches
   - rebuild docs after each batch and report major warning changes

## Review Standard

Every page should answer:

- Who is this for?
- What need does it serve?
- Which Diataxis mode is it in?
- Is that mode pure, or is the page mixed?
- Can a user find it quickly from the homepage/sidebar?

If a page is mixed, prefer split-or-relocate over stuffing more content into it.
