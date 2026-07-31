"""Pydantic models for LabCraft experiments.

Modèles Pydantic pour les expériences LabCraft.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from labcraft.schemas.validators import is_valid_dna_iupac


Chemistry = Literal["LAMP", "RT-LAMP", "PCR", "RT-PCR"]
Enzyme = Literal["Bst", "Bst2.0", "Bst2.0_WarmStart", "Taq"]
ChelationModel = Literal["simple", "full"]
DuplexBackendType = Literal["native", "primer3", "vienna"]
TargetType = Literal["DNA", "RNA"]


class ExperimentConfig(BaseModel):
    """Experiment configuration.
    Configuration de l'expérience.
    """
    model_config = ConfigDict(strict=False)

    name: str
    chemistry: Chemistry
    temperature_C: float = Field(gt=0.0, le=100.0)
    enzyme: Enzyme


class BufferConfig(BaseModel):
    """Buffer configuration.
    Configuration du tampon.
    """
    model_config = ConfigDict(strict=False)

    Na_mM: float = Field(default=50.0, ge=0.0)
    K_mM: float = Field(default=0.0, ge=0.0)
    Tris_mM: float = Field(default=0.0, ge=0.0)
    Mg_total_mM: float = Field(ge=0.0)
    dNTP_total_mM: float = Field(ge=0.0)
    EDTA_mM: float = Field(default=0.0, ge=0.0)
    chelation_model: ChelationModel = Field(default="simple")


class TargetConfig(BaseModel):
    """Target configuration.
    Configuration de la cible.
    """
    model_config = ConfigDict(strict=False)

    id: str
    sequence_file: str
    type: TargetType = Field(default="DNA")
    copies_per_uL: float = Field(default=1000.0, gt=0.0)


class ExplicitDomains(BaseModel):
    """Explicit domain definitions for primers.
    Définitions de domaines explicites pour les amorces.
    """
    model_config = ConfigDict(strict=False)

    F1c: str | None = None
    F2: str | None = None
    B1c: str | None = None
    B2: str | None = None
    linker: str | None = None


DomainSpec = Literal["auto"] | ExplicitDomains


class PrimerConfig(BaseModel):
    """Primer configuration.
    Configuration de l'amorce.
    """
    model_config = ConfigDict(strict=False)

    seq: str
    conc_uM: float = Field(gt=0.0)
    domains: DomainSpec | None = None

    @field_validator("seq")
    @classmethod
    def validate_seq(cls, v: str) -> str:
        """Validate sequence against DNA IUPAC characters.
        Valide la séquence avec les caractères ADN IUPAC.
        """
        if not is_valid_dna_iupac(v):
            raise ValueError(f"Sequence '{v}' contains invalid DNA IUPAC characters.")
        return v


class PrimerSetConfig(BaseModel):
    """Primer set configuration.
    Configuration du jeu d'amorces.
    """
    model_config = ConfigDict(strict=False)

    target: str
    primers: dict[str, PrimerConfig]


class AnalysisConfig(BaseModel):
    """Analysis configuration.
    Configuration de l'analyse.
    """
    model_config = ConfigDict(strict=False)

    temperature_scan: tuple[float, float, float] | None = None
    specificity_db: str | None = None
    amplifiable_dimer_threshold_kcal: float = Field(default=-2.0)
    duplex_backend: DuplexBackendType = Field(default="primer3")
    max_complex_size: int = Field(default=2, ge=1, le=4)


class LabCraftInput(BaseModel):
    """LabCraft main input schema.
    Schéma d'entrée principal LabCraft.
    """
    model_config = ConfigDict(strict=False)

    experiment: ExperimentConfig
    buffer: BufferConfig
    targets: list[TargetConfig]
    primer_sets: list[PrimerSetConfig]
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)

    @model_validator(mode="after")
    def validate_targets_referenced(self) -> LabCraftInput:
        """Validate that targets in primer sets exist.
        Valide que les cibles (targets) des primer sets existent.
        """
        target_ids = {t.id for t in self.targets}
        for pset in self.primer_sets:
            if pset.target not in target_ids:
                raise ValueError(
                    f"Target '{pset.target}' referenced by primer set not found in targets list. "
                    f"Available targets: {list(target_ids)}"
                )
        return self

    @model_validator(mode="after")
    def validate_chemistry_enzyme(self) -> LabCraftInput:
        """Validate chemistry vs enzyme.
        Valide la chimie par rapport à l'enzyme.
        """
        chem = self.experiment.chemistry
        enz = self.experiment.enzyme
        if chem in ("LAMP", "RT-LAMP"):
            if enz not in ("Bst", "Bst2.0", "Bst2.0_WarmStart"):
                raise ValueError(f"For chemistry {chem}, enzyme must be a Bst variant, got {enz}.")
        elif chem in ("PCR", "RT-PCR"):
            if enz != "Taq":
                raise ValueError(f"For chemistry {chem}, enzyme must be Taq, got {enz}.")
        return self

    @model_validator(mode="after")
    def check_mg_dntp(self) -> LabCraftInput:
        """Check Mg >= dNTP (warning).
        Vérifie Mg >= dNTP (avertissement).
        """
        mg = self.buffer.Mg_total_mM
        dntp = self.buffer.dNTP_total_mM
        if mg <= dntp:
            warnings.warn("Mg_total_mM is less than or equal to dNTP_total_mM.")
        return self

    @model_validator(mode="after")
    def validate_primers_for_chemistry(self) -> LabCraftInput:
        """Validate primers for specific chemistries.
        Valide les amorces pour des chimies spécifiques.
        """
        chem = self.experiment.chemistry
        for pset in self.primer_sets:
            names = set(pset.primers.keys())
            if chem in ("LAMP", "RT-LAMP"):
                req = {"F3", "B3", "FIP", "BIP"}
                if not req.issubset(names):
                    raise ValueError(f"For {chem}, primer set must contain at least {req}.")
                
                fip = pset.primers.get("FIP")
                bip = pset.primers.get("BIP")
                if fip and fip.domains is None:
                    raise ValueError("For LAMP, FIP must have domains defined (auto or explicit).")
                if bip and bip.domains is None:
                    raise ValueError("For LAMP, BIP must have domains defined (auto or explicit).")
                
                valid_lamp = {"F3", "B3", "FIP", "BIP", "LF", "LB"}
                invalid = names - valid_lamp
                if invalid:
                    raise ValueError(f"Invalid primer names for {chem}: {invalid}")
            elif chem in ("PCR", "RT-PCR"):
                valid_pcr = {"Forward", "Reverse", "Probe"}
                invalid = names - valid_pcr
                if invalid:
                    raise ValueError(f"Invalid primer names for {chem}: {invalid}")
        return self
        
    @classmethod
    def from_yaml(cls, path: str | Path) -> LabCraftInput:
        """Load and validate from a YAML file.
        Charge et valide depuis un fichier YAML.
        """
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return cls.model_validate(data)
