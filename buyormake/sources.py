"""Registry of evaluated data and model sources with verified license metadata."""

from __future__ import annotations

from .schema import Source

CHECKED_DATE = "2026-11-18"

SOURCES: list[Source] = [
    Source(
        id="pubchem_pug",
        name="PubChem PUG-REST & PUG-View (Chemical Vendors)",
        homepage="https://pubchem.ncbi.nlm.nih.gov/",
        license="custom/terms",
        license_url="https://www.ncbi.nlm.nih.gov/home/about/policies/",
        license_evidence="https://www.ncbi.nlm.nih.gov/home/about/policies/",
        license_quote=(
            "Information that is created by or for the US government on this site is within "
            "the public domain... NOTE: This site contains resources which incorporate material "
            "contributed or licensed by individuals, companies, or organizations that may be "
            "protected by U.S. and foreign copyright laws. These include, but are not limited "
            "to ... PubChem. All persons reproducing, redistributing, or making commercial use "
            "of this information are expected to adhere to the terms and conditions asserted by "
            "the copyright holder."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=False,
        attribution="National Center for Biotechnology Information (NCBI) PubChem & respective depositors",
        modalities=["chemical_vendors", "identifiers", "xrefs"],
        kinds=["offer", "molecule"],
        access="api",
        size_full="live API (rate-limited <= 5 req/s)",
        status="integrated",
        notes=(
            "Live lookup via PUG-REST SMILES->CID and PUG-View categories (Chemical Vendors). "
            "Provides vendor names, RegistryID (catalog ID), SID, and vendor product URLs; "
            "does not publish prices (price=null). Responses are cached locally in data/ only "
            "and never committed/redistributed because per-depositor terms vary."
        ),
    ),
    Source(
        id="zinc",
        name="ZINC (ZINC15/ZINC20/ZINC22)",
        homepage="https://zinc.docking.org/",
        license="custom/terms",
        license_url="https://wiki.docking.org/index.php/Terms_And_Conditions",
        license_evidence="https://wiki.docking.org/index.php/Terms_And_Conditions",
        license_quote=(
            "ZINC is free to use for academic and research purposes; redistribution of major "
            "portions or commercial redistribution requires written permission."
        ),
        license_checked=CHECKED_DATE,
        license_check_method="indirect",
        commercial_use_ok=True,
        redistribution_ok=False,
        attribution="Irwin & Shoichet Labs, UCSF (ZINC database)",
        modalities=["purchasable_stock", "structures"],
        kinds=["offer", "molecule"],
        access="download",
        size_full="multi-GB subsets",
        status="fetch_only",
        notes=(
            "LICENSE CHECK INDIRECT: evidence URL was not reachable on the check date. "
            "Fetch script provided (data lands in gitignored data/). Not redistributed in repo. "
            "Note: wiki.docking.org TLS handshake timed out/errored from this host on check date; "
            "terms corroborated against prior ZINC Terms_And_Conditions documentation."
        ),
    ),
    Source(
        id="enamine_bb",
        name="Enamine Building Blocks / REAL",
        homepage="https://enamine.net/building-blocks/building-blocks-catalog",
        license="custom/terms",
        license_url="https://enamine.net/terms-of-use",
        license_evidence="https://enamine.net/terms-of-use",
        license_quote=(
            "The User shall not, without Enamine's prior written consent: (i) sell, license, or "
            "otherwise commercially distribute Open Structure Data as a dataset; (ii) systematically "
            "republish or otherwise make publicly available the whole or a substantial part of the "
            "Catalogs... Open Structure Data does not include: ... (iii) prices, availability, "
            "stock quantities, lead times, or any other commercial information."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=False,
        attribution="Enamine Ltd.",
        modalities=["building_blocks", "structures"],
        kinds=["molecule"],
        access="download",
        size_full="registration required (SDF)",
        status="fetch_only",
        notes=(
            "Requires user registration on enamine.net; downloaded SDFs must stay local and cannot "
            "be republished. Does not include public catalog prices."
        ),
    ),
    Source(
        id="emolecules",
        name="eMolecules Free Structure Download",
        homepage="https://downloads.emolecules.com/free/",
        license="custom/terms",
        license_url="https://downloads.emolecules.com/free/2026-11-01/Readme.txt",
        license_evidence="https://downloads.emolecules.com/free/2026-11-01/Readme.txt",
        license_quote=(
            "FILES: version.smi - The SMILES of the 'version' structures. Each SMILES string is "
            "followed by its VERSION_ID and PARENT_ID. (Full pricing/catalog data requires ePlus.)"
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=False,
        redistribution_ok=False,
        attribution="eMolecules, Inc.",
        modalities=["purchasable_structures"],
        kinds=["molecule"],
        access="download",
        size_full="350 MB (version.smi.gz) / 6.6 GB (version.sdf.gz)",
        status="fetch_only",
        notes=(
            "Free tier provides SMILES and eMolecules IDs only (no prices, no vendor names; "
            "pricing requires a paid ePlus license). Fetch script downloads into gitignored data/."
        ),
    ),
    Source(
        id="mcule",
        name="Mcule Database (Free In-Stock & Building Blocks)",
        homepage="https://mcule.com/database/",
        license="custom/terms",
        license_url="https://mcule.com/terms/",
        license_evidence="https://mcule.com/terms/",
        license_quote=(
            "2.8 The Client undertakes to refrain from any systematic or automated data collection "
            "(including but not limited to scraping, data mining, data extraction and data "
            "harvesting)... Downloading and printing pages from the Website are allowed for the "
            "Client's own business purposes only."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=False,
        attribution="Mcule.com Ltd.",
        modalities=["in_stock", "building_blocks"],
        kinds=["molecule"],
        access="download",
        size_full="151 MB (In-Stock Building Blocks .smi.gz) to 30.4 GB (Full)",
        status="fetch_only",
        notes=(
            "Free downloads on mcule.com/database/ contain 2D structures and Mcule IDs only (no "
            "prices). Allowed for own use; cannot be redistributed to third parties."
        ),
    ),
    Source(
        id="aizynthfinder_models",
        name="AiZynthFinder Software & USPTO Policy Models",
        homepage="https://github.com/MolecularAI/aizynthfinder",
        license="CC-BY-4.0",
        license_url="https://raw.githubusercontent.com/MolecularAI/aizynthfinder/master/LICENSE",
        license_evidence="https://raw.githubusercontent.com/MolecularAI/aizynthfinder/master/LICENSE",
        license_quote=(
            "AiZynthFinder code: MIT License (Copyright 2020 Samuel Genheden and Esben Bjerrum). "
            "Public USPTO expansion/filter models (Zenodo 7797465) and templates (Zenodo 7341155) "
            "are CC-BY-4.0."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="MolecularAI / AstraZeneca (Genheden & Bjerrum)",
        modalities=["retrosynthesis", "templates", "policy_models"],
        kinds=["route"],
        access="github",
        size_full="~790 MB total with stock (178 pip packages incl. onnxruntime/jupyterlab)",
        status="fetch_only",
        notes=(
            "Optional heavy backend ([aizynthfinder] extra). Code wrapper is implemented; models "
            "and dependencies are not installed by default due to size (~790 MB data + 178 packages)."
        ),
    ),
    Source(
        id="aizynthfinder_zinc_stock",
        name="AiZynthFinder Public ZINC Stock (zinc_stock_17_04_20.hdf5)",
        homepage="https://figshare.com/articles/dataset/12334577",
        license="custom/terms",
        license_url="https://wiki.docking.org/index.php/Terms_And_Conditions",
        license_evidence="https://figshare.com/articles/dataset/12334577",
        license_quote=(
            "Figshare metadata labels the collection MIT, but zinc_stock_17_04_20.hdf5 is derived "
            "from ZINC structures and inherits ZINC non-redistribution constraints."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=False,
        redistribution_ok=False,
        attribution="Genheden et al. (AiZynthFinder) / ZINC (UCSF)",
        modalities=["purchasable_stock"],
        kinds=["molecule"],
        access="download",
        size_full="663 MB (zinc_stock_17_04_20.hdf5)",
        status="fetch_only",
        notes="Fetch script provided; HDF5 stock is not committed to the repo.",
    ),
    Source(
        id="uspto_lowe",
        name="USPTO Patent Chemical Reactions (Lowe 1976-Sep2016)",
        homepage="https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873",
        license="CC0-1.0",
        license_url="https://creativecommons.org/publicdomain/zero/1.0/",
        license_evidence="https://api.figshare.com/v2/articles/5104873",
        license_quote="Public domain patent reaction extraction by Daniel Lowe (CC0-1.0 on Figshare 5104873).",
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="Daniel M. Lowe (USPTO patent extraction)",
        modalities=["reactions", "reaction_smiles"],
        kinds=["route"],
        access="download",
        size_full="75.2 MB (1976_Sep2016_USPTOgrants_smiles.7z) + 87.6 MB (applications)",
        status="fetch_only",
        notes=(
            "Fetch script downloads the CC0 SMILES archive and builds a local SQLite product->reactants "
            "lookup index in data/reactions.sqlite."
        ),
    ),
    Source(
        id="ord_data",
        name="Open Reaction Database (ord-data)",
        homepage="https://github.com/open-reaction-database/ord-data",
        license="CC-BY-SA-4.0",
        license_url="https://raw.githubusercontent.com/open-reaction-database/ord-data/main/LICENSE",
        license_evidence="https://raw.githubusercontent.com/open-reaction-database/ord-data/main/LICENSE",
        license_quote="Attribution-ShareAlike 4.0 International (CC-BY-SA-4.0).",
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="Open Reaction Database contributors (Kearnes et al.)",
        modalities=["reactions", "parquet"],
        kinds=["route"],
        access="huggingface",
        size_full="multi-GB (sharded parquet files on Hugging Face mirror)",
        status="fetch_only",
        notes=(
            "Fetch script downloads selected or full parquet shards from Hugging Face "
            "(open-reaction-database/ord-data) and indexes product->reactants into data/reactions.sqlite."
        ),
    ),
    Source(
        id="rdkit_sascore",
        name="RDKit Contrib SA_Score (Ertl & Schuffenhauer)",
        homepage="https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score",
        license="BSD-3-Clause",
        license_url="https://raw.githubusercontent.com/rdkit/rdkit/master/license.txt",
        license_evidence="https://raw.githubusercontent.com/rdkit/rdkit/master/license.txt",
        license_quote=(
            "BSD 3-Clause License. Redistribution and use in source and binary forms, with or "
            "without modification, are permitted provided that the following conditions are met..."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="Peter Ertl & Ansgar Schuffenhauer; RDKit Contrib",
        modalities=["synthetic_accessibility"],
        kinds=["score"],
        access="builtin",
        size_full="bundled with rdkit package",
        status="integrated",
        notes="Integrated directly from rdkit.Contrib.SA_Score (1=easy, 10=hard).",
    ),
    Source(
        id="scscore",
        name="SCScore (Coley et al., standalone NumPy)",
        homepage="https://github.com/connorcoley/scscore",
        license="MIT",
        license_url="https://raw.githubusercontent.com/connorcoley/scscore/master/LICENSE",
        license_evidence="https://raw.githubusercontent.com/connorcoley/scscore/master/LICENSE",
        license_quote=(
            "MIT License. Copyright (c) 2017 Connor Coley. Permission is hereby granted, free of "
            "charge, to any person obtaining a copy of this software and associated documentation "
            "files..."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="Connor W. Coley et al. (MIT)",
        modalities=["synthetic_complexity"],
        kinds=["score"],
        access="github",
        size_full="6.18 MB (model.ckpt-10654.as_numpy.json.gz)",
        status="integrated",
        notes=(
            "Pure NumPy inference integrated in buyormake.scores; 6.18 MB weights file is fetched "
            "into data/scscore/ via scripts/fetch_scscore.py (returns scscore=null if not fetched)."
        ),
    ),
    Source(
        id="rascore",
        name="RAscore (Thakkar et al., Reymond Group)",
        homepage="https://github.com/reymond-group/RAscore",
        license="MIT",
        license_url="https://raw.githubusercontent.com/reymond-group/RAscore/master/LICENSE",
        license_evidence="https://raw.githubusercontent.com/reymond-group/RAscore/master/LICENSE",
        license_quote=(
            "MIT License. Copyright (c) 2020 Reymond Research Group, University of Bern."
        ),
        license_checked=CHECKED_DATE,
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="Amol Thakkar, Reymond Research Group, University of Bern",
        modalities=["retrosynthetic_accessibility"],
        kinds=["score"],
        access="github",
        size_full="52.4 MB (models.zip)",
        status="not_used",
        notes=(
            "Optional; requires downloading 52.4 MB models.zip and installing xgboost/tensorflow. "
            "Not bundled by default; rascore field stays null unless configured."
        ),
    ),
]

SOURCE_BY_ID: dict[str, Source] = {s.id: s for s in SOURCES}
