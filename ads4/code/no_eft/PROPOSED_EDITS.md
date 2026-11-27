# Proposed edits (NOT applied): replacing the bulk-EFT clause of A2 for the cubic results

Each edit gives the file, the exact OLD text (copied from the file as of 27 Nov 2026) and the NEW text. The justification is
`NOTE.md`, by section number. Nothing here has been applied.

Two styles are possible for the replacement assumption:

* **(S1) recommended.** State the functional-independent assumption (ii′_E), "detector decoupling", and its OPE-data form (ii′).
  The light-remainder step is then graded honestly [D-sketch] through the light-state action bound (new row L0).
* **(S2) alternative.** State the functional-weighted form (ii′_Ω) as the assumption. That keeps item (1) at [D given A2′], but the
  assumption then mentions Ω^EE.

The edits below implement (S1). (S2) changes only E2's statement and the grades of E3 and E10.

---

## E1. letter.tex, assumption (ii) (lines 88–90)

OLD:
```
zero and two may be light. (ii) \emph{Weak coupling below the gap}: $c_T\to\infty$ at fixed $\Delta_g$,
and the states below the gap are described by a tree-level bulk effective theory, uniformly in spin. All
bounds hold at leading order in $1/c_T$. (iii) \emph{Positivity}: the functionals introduced below are
```
NEW:
```
zero and two may be light. (ii) \emph{Weak coupling below the gap}: $c_T\to\infty$ at fixed $\Delta_g$, and below the
gap the OPE of the unit-normalized $T$ with itself is that of a generalized free field up to corrections that vanish
uniformly in spin: the double traces of $T$ have anomalous dimensions $o(c_T^{-1/2})$, and every other state below
the gap enters with squared OPE coefficient $o(1/c_T)$. Large-$N$ factorization gives this at each fixed spin; the
uniformity in spin is the assumption, and no bulk effective theory is assumed for the bounds on $\tau$. All
bounds hold at leading order in $1/c_T$. (iii) \emph{Positivity}: the functionals introduced below are
```
(NOTE §3.3, (ii′). The quartic paragraph's "strongest" input, lines 96–98, already names the tree-level bulk theory for g₄. It is
unchanged. See E1b.)

## E1b. letter.tex (lines 96–98), optional precision

OLD:
```
standard inputs are listed and graded in the Supplemental Material~\cite{SM}. The strongest, used only to identify
the coefficient the quartic sum rule measures with $g_4$, is that the light double-trace data are those of a
tree-level bulk theory, with no fine-tuned cancellation between the heavy states and the tower of bulk contact terms.
```
NEW:
```
standard inputs are listed and graded in the Supplemental Material~\cite{SM}. The strongest, used only for the quartic
bound, is that the light double-trace data are those of a
tree-level bulk theory, with no fine-tuned cancellation between the heavy states and the tower of bulk contact terms.
```
(NOTE §2.2. At k = 4 the double traces enter at first order, so even the existence of the k = 4 light side as a limit needs more
than (ii).)

## E1c. letter.tex, "How the bounds arise" (line 195), optional

OLD:
```
At leading order the light side is fixed by the bulk couplings: graviton exchange, $|\hat g_3|^2$ and, for
the quartic functional, $g_4$.
```
NEW:
```
At leading order the light side of the cubic functional is the action of the stress-tensor block alone, fixed by $c_T$
and $\tau$, because the double traces of $T$ sit at double zeros of the functional; for the quartic functional it also
contains $g_4$.
```

---

## E2. supplement.tex, Table `tab:inputs`, row A2 (line 329)

OLD:
```
A2 & $c_T\to\infty$ at fixed $\Delta_g$; light data ($\Delta-J<\Delta_g$, except light single traces) $=$ tree-level EFT $+\,o(1/c_T)$ uniformly in spin; sum rules used at finite $c_T$, heavy side never expanded; cf.\ CHMRS's tree-level suppression $\varepsilon$ ``much smaller than any power of $1/\Delta_{\rm gap}$'' (Ref.~\cite{CHMRS}, Sec.~1) & D, Q & [A]\\
```
NEW:
```
A2$'$ & $c_T\to\infty$ at fixed $\Delta_g$ (discrete spectrum, finitely many light single traces); sum rules used at finite $c_T$, heavy side never expanded; \emph{detector decoupling}: $\mathcal E_{\rm light}:=\sum_{\mathcal O\in\mathcal L'}(1+J)^{5/2}\log(2+J)\sum_b w_{\mathcal O,b}|a_{\mathcal O,b,0}|^2\to0$, with $\mathcal L'$ the states of twist $<\Delta_g$ other than $1$, $T$ and light single traces, $a_{\mathcal O,b,0}$ the detector amplitudes and $w_{\mathcal O,b}=[T_\Delta-b(b-1)]/3+\tfrac32$ the endpoint weights of Sec.~\ref{sec:angular} (unweighted, $\mathcal E_{\rm light}$ is the light share of the endpoint sum rule, $\le3$ at every $c_T$); implied (two-scale crossing, Sec.~\ref{sec:angcross}) by: $c_T\gamma_{n,J}^2\to0$ uniformly on the light $[TT]_{n,J}$ and squared OPE coefficients (unit-normalised $T$) $o(1/c_T)$ for all other light states and structures, uniformly in spin shells. No bulk effective theory; at each fixed spin it is large-$N$ factorisation. Q also needs L7(b) & D, Q & [A] (the uniformity in spin); fixed spin [D]\\
```

## E3. supplement.tex, row H_K (line 338), last cell

OLD:
```
& D, Q & [D]; order of limits [D given A2]\\
```
NEW:
```
& D, Q & [D]; order of limits [D given A2$'$], light remainder [D-sketch] through L0\\
```

## E4. supplement.tex, new row L0 (insert before the row L1, line 339)

NEW (inserted line):
```
L0 & stripped action of $\Omega^{EE}$ on light states ($\Delta-J<\Delta_g$), Eq.~\eqref{eq:fp-target} continued: $\omega[\mathcal O]=\tfrac{2\pi G\Lambda^4}9\sum_b|a_{\mathcal O,b,0}|^2\hat h_b$, $|\hat h_b|\le T_\Delta[C_h+C'(1+J)^{5/2}\log(2+J)]$, $C_h\approx e^{\pi\Lambda/2}\,10^{4.5}$ ($10^{685.6}$ at $\Delta_g=1000$, $10^{209}$ at $300$); no sign & D & [NE] ($J\le1000$, ball arithmetic with the moment balls of Theorem FP); continuation to light states and $J\to\infty$ asymptotics [D-sketch]\\
```

## E5. supplement.tex, row L1 (line 339), first cell

OLD:
```
L1 & $EE$ admissible; $[TT]$ at double zeros, contacts drop ($k=2,3$ unsubtracted);
```
NEW:
```
L1 & $EE$ admissible; $[TT]$ at double zeros, so at $k=2,3$ (unsubtracted) the $O(1/c_T)$ double-trace data, exchange and contact alike, do not contribute, and their $O(c_T^{-2})$ remainder summed over spin is controlled by A2$'$;
```

## E6. supplement.tex, row L3 (line 356), grade cell

OLD:
```
& D, Q & $c_1=0$ [D] given L1's double zeros (both parities) and L4;
```
NEW:
```
& D, Q & $c_1=0$ [D] given L1's double zeros (both parities), with the light side evaluated block-wise ($T$ block, Lemma of Sec.~\ref{sec:HK}); L4 only as the identity relating the $T$-block action to the Regge arc of the exchange diagram;
```

## E7. supplement.tex, row L4 (line 357), `used` column

OLD:
```
at $n_X=k$ a linear functional of the top Regge coefficient & Q &
```
NEW:
```
at $n_X=k$ a linear functional of the top Regge coefficient & Q (D only as an identity on the exchange diagram, cf.\ L3) &
```

## E8. supplement.tex, row L7 (line 361)

OLD:
```
(b) light double-trace data $=$ tree-level EFT, no conspiracy in heavy sum / contact tower (cf.\ the technical assumption, eq.~(5.3) of Ref.~\cite{CHMRS}) & Q &
```
NEW:
```
(b) light double-trace data $=$ tree-level EFT, no conspiracy in heavy sum / contact tower (cf.\ the technical assumption, eq.~(5.3) of Ref.~\cite{CHMRS}); used for the existence of the $O(1/c_T)$ light side at $k=4,5$, where double traces have single zeros (A2$'$ does not give this), and for $\tilde g_4=g_4$ & Q &
```

## E9. supplement.tex, "Order of limits" (lines 677–686)

OLD:
```
\emph{Order of limits [D given A2].} $H_K$ is used only for the exact correlator at finite $c_T$, where R1 holds. (Pos) is a statement
about $(\Delta,J)$ and parity, so it applies at every $c_T$ to every state with $\Delta-J\ge\Delta_g$, single- or multi-trace, and to
light single traces; the sum rule then gives $\Omega[{\rm light}](c_T)\le0$. At $k=2,3$ the light side is $c_T^{-1}$ times the action
on the $T$ block with the exact $\langle TTT\rangle$ coefficients, plus light double and multi traces, which enter at $O(c_T^{-2})$
(double zeros). At $k=4,5$ the first double-trace families enter at $O(c_T^{-1})$ with their tree-level data (L4--L7). Given A2,
$\Omega[{\rm light}]=c_T^{-1}L_1+o(c_T^{-1})$, so $L_1\le0$, which is the Result.
```
NEW:
```
\emph{Order of limits [D given A2$'$; light remainder D-sketch through L0].} $H_K$ is used only for the exact correlator at finite $c_T$, where R1 holds. (Pos) is a statement
about $(\Delta,J)$ and parity, so it applies at every $c_T$ to every state with $\Delta-J\ge\Delta_g$, single- or multi-trace, and to
light single traces; the sum rule then gives $\Omega[{\rm light}](c_T)\le0$. At $k=2,3$ the light side is exactly
$c_T^{-1}L_1+R(c_T)$: $c_T^{-1}L_1$ is the action of the single $T$ block with the exact $(c_T,t_4,\alpha_T)$, a conformal-kinematic
quantity with no $1/c_T$ expansion [D], and $R$ is the sum over the light states $\mathcal L'$ other than $1$, $T$ and light single
traces. Hence, at every $c_T$ and with no light single trace, $\tau(c_T)^2\le\tau_0^2+(57600/K\Lambda^8)\,\pi c_T|R(c_T)|/(32n_0)$ [P],
$n_0>0$ the normalisation of the $T$-block term, and "leading order in $1/c_T$" means $\limsup_{c_T\to\infty}\tau\le\tau_0$, which holds
iff $c_TR\to0$. By L0, $c_T|R|\le(2\pi c_TG\Lambda^4/3)\max(C_h,C')\,\mathcal E_{\rm light}(c_T)$ up to $O(1)$ factors, so A2$'$ gives
$c_TR\to0$; at each fixed spin the double zeros make the contribution $O(c_T^{-2})$ under factorisation alone, and only the uniformity
in spin is assumed. No bulk effective theory and no Witten diagram enters at $k=2,3$. At $k=4,5$ the first double-trace families
enter at $O(c_T^{-1})$ (single zeros), and the existence of the light side and its value need L4--L7, including L7(b).
```
and, in the same paragraph,

OLD:
```
\cite[Sec.~2.2.2]{CLSD} but do not discuss the order of limits. The uniformity in spin of the light-side
remainder is [D-sketch] inside A2, and its size (a power $\Delta_g^p/c_T^2$) is not estimated: the statements hold as $c_T\to\infty$, with no statement on how large $c_T$ must be for the numbers to bound the exact $t_4$, $\alpha_T$, $\tilde g_4$.
```
NEW:
```
\cite[Sec.~2.2.2]{CLSD} but do not discuss the order of limits. The size of the light-side remainder is not a power of $\Delta_g$:
term by term it is $O(C_h/c_T^2)$ times a polynomial in $\Lambda$, with $C_h\approx e^{\pi\Lambda/2}10^{4.5}\approx10^{686}$ at
$\Delta_g=1000$ (L0), so the inequality above bounds the exact $t_4$, $\alpha_T$ only for $c_T\gtrsim10^{690}$, unless the spin sum of the
light double traces cancels; this holds with or without a bulk effective theory. The statements hold as $c_T\to\infty$.
```

## E10. supplement.tex, Table `tab:c1`, rows (i), (ii) (lines 794–795)

OLD:
```
(i) Lemma-R value of a tree term $=$ its term-wise $s$-channel value (L4) & [D] (inherited)\\
(ii) at $O(1/c_T)$ only the $T$ block survives: $L[T]|\Omega\rangle=0$; $[TT]$ at double zeros (L1), both parities & [D] (inherited; parity-odd zeros in closed form)\\
```
NEW:
```
(i) the light side is the block-wise $s$-channel sum (definition, Sec.~\ref{sec:arc}); the $T$-block action equals the Regge arc of the exchange diagram (L4, an identity on that function; not needed if the $T$ block is evaluated directly, Steps 1--3) & [D]\\
(ii) light side $=$ $T$ block $+\,R$, $L[T]|\Omega\rangle=0$; $[TT]$ at double zeros (L1), both parities; $c_TR\to0$ by A2$'$ & [D] (zeros; parity-odd in closed form); $c_TR\to0$ [D-sketch given A2$'$] (L0)\\
```
The caption's "The net grade is that of the inherited links (i), (ii)" can stay. c₁ = 0 itself is a statement about the T-block
term and is exact at every c_T.

## E11. supplement.tex, angular route (lines 2112–2114) and the Letter-mapping paragraph (line 106)

OLD (2112–2114):
```
$\mathcal L_g\to0$ follows from
planar factorisation: transitions from $T$ to multi-trace states are $O(c_T^{-1/2})$, and $\Delta_T$ is protected, so degenerate
multi-trace clusters have no $O(c_T^{-1/2})$ splitting [D].
```
NEW:
```
$\mathcal L_g\to0$ follows from
A2$'$: by the ladder ratio, $|a_{\mathcal O,4,2}|^2\le101|a_{\mathcal O,4,0}|^2$ on light states and $w_{\mathcal O,4}\ge\tfrac32$, so
$\mathcal L_g\le70\,\mathcal E_{\rm light}\to0$ [D given A2$'$]. State by state this is planar factorisation (transitions from $T$ to
multi-trace states are $O(c_T^{-1/2})$, and the double zeros make degenerate clusters harmless). The sum over the infinitely many
light spins is the content of A2$'$.
```
OLD (106):
```
$\mathcal L_g\to0$ [D] from planar factorisation;
```
NEW:
```
$\mathcal L_g\to0$ [D given A2$'$] (planar factorisation state by state, uniformity in spin from A2$'$);
```

## E12. supplement.tex, Remark of Sec. `sec:allgap`.A (lines 1597–1599)

OLD:
```
\emph{Remark} [D]. A theory that satisfies A1$'$ and A2 with gap $\Delta_g'$ satisfies them with every $\Delta_g\le\Delta_g'$:
its single traces of spin $\ge4$ have twist $\ge\Delta_g'\ge\Delta_g$, and the states of twist $<\Delta_g$ are among those of
twist $<\Delta_g'$, whose data are tree level.
```
NEW:
```
\emph{Remark} [P]. A theory that satisfies A1$'$ and A2$'$ with gap $\Delta_g'$ satisfies them with every $\Delta_g\le\Delta_g'$:
its single traces of spin $\ge4$ have twist $\ge\Delta_g'\ge\Delta_g$, and the light set at $\Delta_g$ is a subset of the light set
at $\Delta_g'$, so the sum defining $\mathcal E_{\rm light}$ only loses terms; the states of twist in $[\Delta_g,\Delta_g')$ are heavy
at $\Delta_g$ and are covered by (Pos).
```

## E13. supplement.tex, Net grades after Table `tab:inputs` (line 374) and Sec. `sec:full`

OLD (line 374):
```
\textbf{Net grades} (repeated in Sec.~\ref{sec:full}). Item (1): [D] if no light single trace couples to $TT$ (the Letter's headline);
```
NEW:
```
\textbf{Net grades} (repeated in Sec.~\ref{sec:full}). Item (1): [D] given A1$'$, A2$'$ and the inputs if no light single trace couples to $TT$ (the Letter's headline), with the light-remainder step [D-sketch] through L0 (formerly hidden inside A2);
```
(Repeat the same insertion in the "Net grades" paragraph of Sec. `sec:full`, "Implication (1): [D] if no light single trace
couples to TT", and in the Einstein-rigidity grades. For item (3), add "𝓛_g → 0 [D given A2′]".)

## E14. supplement.tex, the Glossary/Letter-mapping paragraph (line 82)

OLD:
```
\emph{The Letter's assumptions.} (i) and (ii) are A1$'$ and A2 below.
```
NEW:
```
\emph{The Letter's assumptions.} (i) and (ii) are A1$'$ and A2$'$ below; the Letter's wording of (ii) is the sufficient
OPE-data form of A2$'$ (anomalous dimensions $o(c_T^{-1/2})$, other squared OPE coefficients $o(1/c_T)$, uniformly in spin).
```

## E15. Global

Replace "A2" by "A2′" in:
* Sec. `sec:full` (Result, Einstein rigidity);
* Theorem L′ (line 1484): A2 is idle there and could simply be dropped;
* Sec. `sec:cor` (line 1967);
* Sec. `sec:abjm` (line 2187: "satisfies A1′ and A2′ … at leading order, given that the planar light double-twist anomalous
  dimensions are O(Δ_g^p/c_T) uniformly in spin [D-sketch]");
* Scope (line 2202).

Add one line to the Scope: "a bound on the size of the 1/c_T remainder below c_T ≈ 10⁶⁹⁰ (row L0)".

## Not proposed

* No change to rows R1–R3, DS, LRC, L2, L5, L6, L6′, L8, L9, LS1–LS6.
* The side observation of NOTE §6 (light J = 0, 2 block actions from the ladder formula, reproducing F_D⁽⁰⁾ with the mass shift
  δ(Δ_L) = 0) is not proposed for the paper until it is checked against a direct evaluation of the block action.
