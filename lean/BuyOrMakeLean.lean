/-!
# BuyOrMakeLean: Formal L5 Stoichiometric, Strict-Null Budget & Green Pareto Guardrails

Self-contained Lean 4 (`v4.36.0`) formalization for `buy-or-make` (`L5`), synchronized
with the central 14-repository Auto-Lab Trust Regime (`AutoLabLean`):

1. **Acyclic Disconnection & Heavy-Atom Conservation (`AtomVec`, `routeStoichiometryOK`)**:
   - Synthesis routes (`isDirectBuy = false`) must cleave only acyclic bonds
     (`bondIsAcyclic = true`, preventing cyclic lactam/lactone atom duplication) and
     conserve exact heavy-atom counts `sumAtoms reactantAtoms = targetAtoms + leavingDelta`.
2. **Strict-Null Homogeneous `USD/g` Budget Gate (`budgetPolicyOK`, `unpriced_never_passes_budget`)**:
   - Missing catalog prices (`priceCents = none`), incomplete leaf pricing
     (`costComplete = false`), or non-`USD/g` units (`isUsdPerGram = false`) can NEVER
     pass an active numeric budget gate (`budgetCents = some B`).
3. **Transitive License Gate & 3D Green Pareto Dominance (`licensePolicyOK`, `GreenTriplet`)**:
   - Propagates license constraints across target stock, reaction DBs, and leaf synthons,
     and defines decidable 3D dominance over `(costCentsPerGram, leadTimeDays, pmiPermille)`.
-/
namespace BuyOrMakeLean

/-- License posture harmonized across `buy-or-make` and the Auto-Lab Trust Regime. -/
inductive LicenseClass where
  | permissive
  | shareAlike
  | nonCommercial
  | evalOnly
  | customTerms
  deriving DecidableEq, Repr

def LicenseClass.commercialTrainOK : LicenseClass → Bool
  | .permissive => true
  | _           => false

/-- Structured boundary violation reporting the repository, rule identifier, and detail. -/
structure RegimeViolation where
  repo   : String
  rule   : String
  detail : String
  deriving DecidableEq, Repr

def RegimeViolation.format (v : RegimeViolation) : String :=
  s!"[{v.repo} :: {v.rule}] {v.detail}"

/-- Return the first violation in a list of `(holds, violation)` checks. -/
def firstViolation : List (Bool × RegimeViolation) → Option RegimeViolation
  | [] => none
  | (true, _) :: rest => firstViolation rest
  | (false, v) :: _ => some v

@[simp] theorem firstViolation_nil : firstViolation [] = none := rfl

@[simp] theorem firstViolation_cons (b : Bool) (v : RegimeViolation)
    (rest : List (Bool × RegimeViolation)) :
    firstViolation ((b, v) :: rest) = none ↔ b = true ∧ firstViolation rest = none := by
  cases b <;> simp [firstViolation]

theorem firstViolation_none_iff (checks : List (Bool × RegimeViolation)) :
    firstViolation checks = none ↔ ∀ p ∈ checks, p.1 = true := by
  induction checks with
  | nil => simp [firstViolation]
  | cons hd tl ih =>
    rcases hd with ⟨b, v⟩
    cases b with
    | false => simp [firstViolation]
    | true => simp [firstViolation, ih]

/-! ## 1. Stoichiometric Atom Conservation, Strict Pricing & License Gate -/

/-- Exact heavy-atom count vector `(C, N, O, S, Hal)`. -/
structure AtomVec where
  c   : Nat
  n   : Nat
  o   : Nat
  s   : Nat
  hal : Nat
  deriving DecidableEq, Repr

def AtomVec.add (a b : AtomVec) : AtomVec :=
  { c := a.c + b.c, n := a.n + b.n, o := a.o + b.o, s := a.s + b.s, hal := a.hal + b.hal }

def sumAtoms : List AtomVec → AtomVec
  | [] => { c := 0, n := 0, o := 0, s := 0, hal := 0 }
  | x :: xs => AtomVec.add x (sumAtoms xs)

/-- Procurement or retrosynthetic proposal checked by `buy-or-make` at Layer L5. -/
structure MakeOrBuyProposal where
  isDirectBuy      : Bool
  steps            : Nat
  maxSteps         : Nat
  bondIsAcyclic    : Bool
  targetAtoms      : AtomVec
  leavingDelta     : AtomVec
  reactantAtoms    : List AtomVec
  allLeavesInStock : Bool
  costComplete     : Bool
  isUsdPerGram     : Bool
  priceCents       : Option Nat
  budgetCents      : Option Nat
  sourceLicenses   : List LicenseClass
  commercialMode   : Bool
  deriving DecidableEq, Repr

def routeStoichiometryOK (p : MakeOrBuyProposal) : Bool :=
  if p.isDirectBuy then
    p.steps == 0
  else
    decide (1 ≤ p.steps ∧ p.steps ≤ p.maxSteps) &&
    p.bondIsAcyclic &&
    p.allLeavesInStock &&
    decide (sumAtoms p.reactantAtoms = AtomVec.add p.targetAtoms p.leavingDelta)

def budgetPolicyOK (p : MakeOrBuyProposal) : Bool :=
  match p.budgetCents with
  | none => true
  | some B =>
    match p.priceCents with
    | none => false
    | some c => p.costComplete && p.isUsdPerGram && decide (c ≤ B)

def licensePolicyOK (p : MakeOrBuyProposal) : Bool :=
  if p.commercialMode then
    p.sourceLicenses.all LicenseClass.commercialTrainOK
  else
    p.sourceLicenses.all (fun l => l != .evalOnly)

/-- Declarative specification of a valid L5 procurement proposal. -/
structure MakeOrBuyOK (p : MakeOrBuyProposal) : Prop where
  stoichOK  : routeStoichiometryOK p = true
  budgetOK  : budgetPolicyOK p = true
  licenseOK : licensePolicyOK p = true

/-- Decidable guardrail checker for `buy-or-make`. -/
def checkBuyOrMake (p : MakeOrBuyProposal) : Option RegimeViolation :=
  firstViolation
    [ (routeStoichiometryOK p,
        { repo := "buy-or-make", rule := "acyclic_atom_balanced_route",
          detail := "route violates acyclic (!@) bond constraint, max_steps, stock leaves, or heavy-atom conservation" }),
      (budgetPolicyOK p,
        { repo := "buy-or-make", rule := "strict_homogeneous_budget",
          detail := "missing/incomplete price or non-USD/g unit cannot pass numeric budget gate" }),
      (licensePolicyOK p,
        { repo := "buy-or-make", rule := "transitive_license_gate",
          detail := "target, reaction-DB, or leaf source license incompatible with requested mode" }) ]

/-- Soundness and completeness of `checkBuyOrMake` wrt `MakeOrBuyOK`. -/
theorem checkBuyOrMake_ok_iff (p : MakeOrBuyProposal) :
    checkBuyOrMake p = none ↔ MakeOrBuyOK p := by
  unfold checkBuyOrMake
  simp only [firstViolation_cons, firstViolation_nil, and_true]
  constructor
  · rintro ⟨h1, h2, h3⟩; exact ⟨h1, h2, h3⟩
  · rintro ⟨h1, h2, h3⟩; exact ⟨h1, h2, h3⟩

/-- Hard boundary theorem 1: Endocyclic ring cleavage (`bondIsAcyclic = false`) is rejected. -/
theorem ring_cleavage_rejected (p : MakeOrBuyProposal)
    (hmake : p.isDirectBuy = false) (hring : p.bondIsAcyclic = false) :
    checkBuyOrMake p ≠ none := by
  intro hok
  have hs := ((checkBuyOrMake_ok_iff p).mp hok).stoichOK
  simp [routeStoichiometryOK, hmake, hring] at hs

/-- Hard boundary theorem 2: Unpriced or incomplete-cost candidates NEVER pass an active budget. -/
theorem unpriced_never_passes_budget (p : MakeOrBuyProposal) (B : Nat)
    (hbud : p.budgetCents = some B) (hunpriced : p.priceCents = none ∨ p.costComplete = false) :
    checkBuyOrMake p ≠ none := by
  intro hok
  have hb := ((checkBuyOrMake_ok_iff p).mp hok).budgetOK
  rcases hunpriced with hp | hc
  · simp [budgetPolicyOK, hbud, hp] at hb
  · cases hpc : p.priceCents <;> simp [budgetPolicyOK, hbud, hpc, hc] at hb

/-- Hard boundary theorem 3: Heavy-atom imbalance across reactants and target+leaving group is rejected. -/
theorem atom_imbalance_rejected (p : MakeOrBuyProposal)
    (hmake : p.isDirectBuy = false)
    (himbal : sumAtoms p.reactantAtoms ≠ AtomVec.add p.targetAtoms p.leavingDelta) :
    checkBuyOrMake p ≠ none := by
  intro hok
  have hs := ((checkBuyOrMake_ok_iff p).mp hok).stoichOK
  simp [routeStoichiometryOK, hmake, himbal] at hs

/-- Hard boundary theorem 4: Non-`USD/g` offers cannot pass a numeric per-gram budget gate. -/
theorem non_usd_per_gram_rejected_under_budget (p : MakeOrBuyProposal) (B : Nat)
    (hbud : p.budgetCents = some B) (hunit : p.isUsdPerGram = false) :
    checkBuyOrMake p ≠ none := by
  intro hok
  have hb := ((checkBuyOrMake_ok_iff p).mp hok).budgetOK
  cases hpc : p.priceCents <;> simp [budgetPolicyOK, hbud, hpc, hunit] at hb

/-! ## 2. 3D Cost-LeadTime-Greenness (`PMI`) Pareto Dominance -/

/-- Exact fixed-point 3D procurement objective `(costCentsPerGram, leadTimeDays, pmiPermille)`. -/
structure GreenTriplet where
  costCentsPerGram : Nat
  leadTimeDays     : Nat
  pmiPermille      : Nat
  deriving DecidableEq, Repr

/-- `a` Pareto-dominates `b` iff `a` is no worse in all 3 axes and strictly better in at least one. -/
def paretoDominates (a b : GreenTriplet) : Bool :=
  decide (a.costCentsPerGram ≤ b.costCentsPerGram ∧
          a.leadTimeDays ≤ b.leadTimeDays ∧
          a.pmiPermille ≤ b.pmiPermille ∧
          (a.costCentsPerGram < b.costCentsPerGram ∨
           a.leadTimeDays < b.leadTimeDays ∨
           a.pmiPermille < b.pmiPermille))

/-- Strict irreflexivity: no procurement option Pareto-dominates itself. -/
theorem paretoDominates_irrefl (a : GreenTriplet) :
    paretoDominates a a = false := by
  unfold paretoDominates
  apply decide_eq_false
  omega

/-! ## 3. Compile-Time `#guard` Witnesses -/

def demoValidRoute : MakeOrBuyProposal :=
  { isDirectBuy := false
    steps := 1
    maxSteps := 2
    bondIsAcyclic := true
    targetAtoms := { c := 12, n := 3, o := 1, s := 0, hal := 0 }
    leavingDelta := { c := 0, n := 0, o := 0, s := 0, hal := 1 }
    reactantAtoms := [ { c := 4, n := 2, o := 0, s := 0, hal := 1 },
                       { c := 8, n := 1, o := 1, s := 0, hal := 0 } ]
    allLeavesInStock := true
    costComplete := true
    isUsdPerGram := true
    priceCents := some 1450
    budgetCents := some 2500
    sourceLicenses := [.permissive]
    commercialMode := true }

#guard (checkBuyOrMake demoValidRoute).isNone == true

-- Endocyclic ring cleavage (`bondIsAcyclic = false`) -> DENIED
#guard (checkBuyOrMake { demoValidRoute with bondIsAcyclic := false }).isSome == true

-- Strict null price (`priceCents = none`) under active budget -> DENIED
#guard (checkBuyOrMake { demoValidRoute with priceCents := none }).isSome == true

-- Incomplete leaf cost (`costComplete = false`) under active budget -> DENIED
#guard (checkBuyOrMake { demoValidRoute with costComplete := false }).isSome == true

-- Heavy-atom imbalance (`c := 11` instead of `12`) -> DENIED
#guard (checkBuyOrMake
         { demoValidRoute with targetAtoms := { c := 11, n := 3, o := 1, s := 0, hal := 0 } }).isSome == true

-- Green catalytic route dominates slow high-waste route
#guard paretoDominates
         { costCentsPerGram := 1800, leadTimeDays := 3, pmiPermille := 14000 }
         { costCentsPerGram := 1800, leadTimeDays := 9, pmiPermille := 65000 }
       == true

end BuyOrMakeLean
