# Mathematical Protocol, Revision 4: Noise/Bandwidth Equivalence and Reference-Guided Correction under Recursive KDE Training

Status: pre-registration draft, revision 4 (supersedes revisions 1–3; earlier content is retained except where Section 0 states a change). Pilot parameter values are carried over unchanged from revision 1 (Section 1) and are **not** final. Nothing marked [U] has been run.

---

## 0. Evidence legend and revision summary

| Tag | Meaning |
|---|---|
| **[A]** | Analytical identity or proof sketch given in the text (not itself a numerical result). |
| **[G]** | Deterministic numerical grid check: population operator iterated on a grid ($dx=0.02$, $x\in[-100,110]$), or exact integration. From `validate.py`, output `validate_out.txt`. |
| **[MC]** | Monte Carlo observation: finite-sample simulation with the stated seeds (same script/output). |
| **[H]** | Hand arithmetic from a formula; not executed in code. |
| **[U]** | Proposal or unvalidated definition; no numerical check exists yet. |

All [G] and [MC] numbers below are copied from `validate_out.txt` without alteration. Items that the validation run did **not** cover are listed in Section 15.

**What changed in revision 2**
1. Rare-mode endpoints are separated into lineage survival/mass, conditional spatial displacement of component-3 descendants, and spatial rare-region mass with a deterministic spill-over comparator (Section 6). Spatial mass is no longer called "mode mass".
2. H1 is reframed as a noise-versus-bandwidth equivalence/negative-control experiment (Section 9). The "noise helps" hypothesis is retired.
3. The horizon is rewritten (Section 7): analytical 5% residual horizon retained, endpoint-specific stabilisation checks added, generation-matched and residual-matched results both reported.
4. Proposition 3 is rewritten in terms of lineage fraction, separating reference-estimation error from recursive deviation (Section 4). The revision-1 statement that extinction is "negligible at $\lambda=0.1$" is **retracted**.
5. A reference-budget factorial (vary $m$ at fixed $\lambda$; vary $\lambda$ at fixed $m$; proportional condition) is added (Section 8).
6. The full validation record is preserved in Section 5 with evidence labels.
7. Hypotheses and endpoints are redefined (Section 12): H1 negative control, H2 primary, AB secondary interaction.

**What changed in revision 4** (response to the revision-3 review; new numerics V11–V14 in Section 5)
1. **$W_1$ promoted to co-primary for H2** (P-H2c; Section 12). H2 is decided on all three of P-H2a, P-H2b, P-H2c with Holm over three; item 23 closed. Numerics give the reference-estimation benchmark for $W_1$ (V13) and show that the pilot spatial error is an *excess* from inflow (V12, V14), which refines the interpretation of P-H2a (Section 6.4).
2. **Seed-conditional population operator** $\mathcal T^{\rm ref}$ added alongside the target-averaged $\mathcal T$ (Section 3); lineage-resolved comparator now has both versions (Section 6.5); $W_1(Q^*_{\rm ref},Q^*)\le W_1(\hat P_{\rm ref},P_0)$ proved and checked (Prop 1, V13).
3. **Proposition 4 withdrawn** (Section 4, 10): its limiting condition was inconsistent with the noise schedule. Replaced by empirical conjectures C1–C3 that carry no weight in any decision.
4. **Interpretation of R for arm B corrected** (Section 7.4): R matches the *uncorrected* operator's residual contraction; it does not guarantee equal convergence stage for B.
5. **Reference draw rule standardised** (Sections 3, 8): one with-replacement rule across F1, F2, F3; the revision-3 balanced-replication/subsampling rule becomes a sensitivity analysis. Exact stationary lineage variance for both rules derived and checked (Prop 3, V11).
6. **Stabilisation rule made sustained** (Section 7.3), with edge cases for $E_*=E_0$, no crossing and oscillation.
7. **"Gap of 6" wording corrected** (Section 9): leakage is cross-region diffusion over ancestry depth with an exact formula, checked against the grid (V14).

**Addendum to revision 4 (validation pass 2; V15–V17; scripts `validate_rev5_*.py`)**
- **Horizon (V15):** $T_5$ is a common exposure horizon, not a stabilisation horizon for the spatial endpoints; $t_{\rm stab}/T_5$ ranges from 0.76 to 2.0 across $(\lambda,h)$. An amended R-horizon rule is **proposed, not applied** (Section 7.4, item 33); the primary R window stays at $T_5$ until you ratify.
- **Full KDE lineage chain (V16):** labels carried through KDE sampling, reference insertion and fitting reproduce $V_{\rm fixed}$ and $V_{\rm wr}$ to within 0.8% at $\lambda=0.10,0.20$ (20,000 chains). The test is an implementation check: label dynamics cannot depend on positions.
- **Transition-matrix controls (V17):** non-transporting arm consistent with $\Delta A=0$ (max $\lvert z\rvert=1.98$ over 9 cells); artificial transporters reproduce the predicted $\Delta A$ (all cells within $\lvert z\rvert\le1.5$) at three strengths, down to a $\Delta A[1,R_3]$ of 17% of the target rare-region mass.
- **$W_1$ interpretation corrected:** $\mathbb E\,W_1(\hat P_{\rm ref},P_0)$ is a benchmark for the size of the reference-estimation discrepancy and the rationale for pairing, **not** a universal error floor (Section 5 V13, Section 14, item 31).
- **Grid resolution:** grid values of $e_R$ carry about $\pm0.003$ from the cell containing $b_{23}$ (V15); closed forms (V14) are the reference.

**What changed in revision 3** (researcher decisions recorded)
1. **Primary horizons designated (Section 7.4):** residual-matched (R) for H2, AB and the factorial F; generation-matched (G) for H1. The other horizon is always reported as secondary and never pooled.
2. **Equivalence margins:** the universal 10% placeholder is retired. Margins are endpoint-specific, each with a stated anchor (Section 9.3). Consequential change: H1 is tested on signed endpoints with a paired-by-reference TOST, so reference-estimation error cancels (Section 9).
3. **Evidence of mass transport:** an ancestry-to-spatial-region transition matrix plus conditional displacement (Section 6.5). It is a **secondary mechanistic endpoint** and has no role in the H2 decision rule: H2 must win on the primary distributional endpoints (Section 12); evidence that samples moved is not evidence of recovery.
4. **Empirical stationarity:** separate, endpoint-specific stabilisation checks with per-endpoint flags and a sensitivity analysis (Section 7.3).
5. Section 15 updated: items 4, 7, 10 resolved or narrowed; new items 21–26 added.

---

## 1. Pilot configuration (carried over from revision 1; not final)

| Item | Pilot value | Note |
|---|---|---|
| Target | $\pi=(0.70,0.27,0.03)$, $\mu=(-3,2,8)$, $\sigma=(1.0,0.8,0.6)$ | |
| Tail region | $R_\tau=(-\infty,-5.5]\cup[9.5,\infty)$ | |
| Training size $N$ | 2000 | |
| Generation count $T$ | 20 (pilot count) | Not an equilibrium horizon; see Section 7 |
| Generation-matched window (G) | generations 16–20 | primary for H1 only |
| Residual-matched window (R) | $[T_5(\lambda)-4,\,T_5(\lambda)]$, Section 7.1 | primary for H2, AB, F |
| Reference fraction $\lambda$ | pilot grid $\{0.01,0.02,0.05\}$ | 0.10 and 0.20 appear in validation tables only |
| Reference size | $m=\lambda N$ (proportional design, all $m$ points used each generation) | the only reference-sampling rule validated |
| Base bandwidth $h_0$ | pilot grid $\{0.3,0.6,0.9\}$ | |
| Noise grid | $\sigma=\kappa h_0/\sqrt{1-\lambda}$, $\kappa\in\{0.3,0.6,0.9\}$, $h(\sigma)^2=h_0^2(1-\kappa^2)$ | |
| Archive size $M_A$ | $10^4$ | |
| Score levels | $s_{\max}\in\{0.5,1,2.5,6\}$; $s_{\min}$ by LOO-CV | |
| Stratified retention | $K=3$ strata, $\gamma\in\{1,0.5,0\}$ | |
| Density-ratio baseline | features $1,x,\dots,x^4$ (or small MLP); $M_s=10^4$; 5-fold cross-fit; $r_{\max}$ pilot | |
| Verifier | $\alpha=0.01$ LOO quantile | |
| Seeds | pilot 1–100; test 1001–1300 ($n=300$) | |
| Tests | one-sided $\alpha=0.025$ per hypothesis; bootstrap $B=10^4$ | |
| Guard margin | 0.05 | |
| Equivalence margins | endpoint-specific, Section 9.3 | replaces the universal 10% placeholder |
| Legacy extinction threshold | $m_R<\tfrac12 m_R^0$ (region mass) | retained only for the pilot observations in V9 |

New in revision 2 (also pilot-only, not run): factorial reference-budget grid, Section 8.

---

## 2. Target distribution, regions, boundaries

$$P_0(x)=\sum_{k=1}^{3}\pi_k\,\varphi(x;\mu_k,\sigma_k^2).$$

**[G]** Mean $\mu_0=-1.32000$; variance $v_0=8.44120$.

**Adjacent-crossing boundaries [G].** $b_{12}$ solves $\pi_1\varphi_1=\pi_2\varphi_2$ on $(-3,2)$; $b_{23}$ solves $\pi_2\varphi_2=\pi_3\varphi_3$ on $(2,8)$:
$$b_{12}=-0.10488,\qquad b_{23}=5.58249.$$
**Spatial regions (protocol definition):** $R_1=(-\infty,b_{12})$, $R_2=[b_{12},b_{23})$, $R_3=[b_{23},\infty)$. Spatial region mass $m_R(Q)=Q(R_3)$.

**Literal argmax issue [G].** The regions $\{x:\ k=\arg\max_j\pi_j\varphi_j(x)\}$ are **not** three intervals. Computed on a $10^{-4}$ grid over $[-80,120]$:

| Interval | argmax component |
|---|---|
| $(-\infty,-0.1049)$ | 1 |
| $[-0.1049,5.5824)$ | 2 |
| $[5.5824,24.3550)$ | 3 |
| $[24.3550,\infty)$ | **1** (largest-variance component dominates the far right tail) |

The grid-resolution difference between 5.5824 and $b_{23}=5.58249$ is the $10^{-4}$ grid step. Under $P_0$ both definitions give the same masses to displayed precision: argmax $[0.69982,0.27018,0.03000]$, adjacent $[0.69982,0.27018,0.03000]$ (component weights $\pi=[0.70,0.27,0.03]$). At the population fixed points in V4, the two definitions differ in $e_R$ by at most 0.0049 (1.7658 vs 1.7609, $\lambda=0.02$ unmatched arm). **The adjacent-crossing definition is the protocol definition**; argmax is reported in V4 for transparency only.

**Tail mass [G].** $\tau(Q)=Q(R_\tau)$; $\tau_0=\tau(P_0)=0.004347+0.000186=0.004533$ (left $0.004347$, right $0.000186$). The tail region is dominated (96%) by the left tail of component 1.

---

## 3. Recursive training operator

**Learner.** $\Phi_h(D,w)=\sum_i w_i\,\varphi(\cdot;x_i,h^2)$ (Gaussian KDE, fixed $h$).

**Labels (new in revision 2) [U].** Every datum $x$ carries a component label $c(x)\in\{1,2,3\}$: archive and $S_0$ points receive the label of the mixture component that generated them (draw $k\sim\pi$, then $x\sim\mathcal N(\mu_k,\sigma_k^2)$; the draw of $k$ must be stored). A KDE draw $y=x_J+h\varepsilon\,(+\sigma\xi)$ inherits $c(y)=c(x_J)$. Selection arms (DR, V) preserve labels. For correction arm B the output point inherits the label of its input point, so the label records **origin**, not position. Because index draws ignore position, the label dynamics of B have the same law as those of C (Prop 3); lineage mass therefore cannot detect transport (Section 6.5).

**Archive and reference.** Per seed, $\mathcal A=\{a_i\}_{i\le M_A}\overset{iid}{\sim}P_0$; retention rule $\rho$ selects $D_{\rm ref}$ ($|D_{\rm ref}|=m$) with weights $\omega_j$ (uniform: $1/m$).

**Reference draw rule.** Each generation uses $n_{\rm ref}=\lambda N$ reference points:
- **Validated design:** $n_{\rm ref}=m$; every reference point is used every generation ([G]/[MC] results in this document all use this rule).
- **With-replacement rule (primary for the factorial, Section 8) [pipeline U; label chain checked, V11]:** every generation draws $n_{\rm ref}=\lambda N$ points **with replacement** from the retained set, for every $m$. One rule across all $m$ and $\lambda$ removes the revision-3 switch between balanced replication ($m<n_{\rm ref}$) and subsampling ($m>n_{\rm ref}$).
- **Revision-3 generalised rule (balanced replication / subsampling without replacement) [U]:** retained only as a secondary sensitivity analysis.

**Generation 0.** $S_0\sim P_0^{\otimes(N-n_{\rm ref})}$ fresh; $D_0=S_0\cup(\text{reference draw})$; $\hat P_0=\Phi_h(D_0)$.
**Generation $t\ge1$.**
1. $Y_t=(y_i)_{i\le N-n_{\rm ref}}\overset{iid}{\sim}\hat P_{t-1}$ (draw $y=x_J+h\varepsilon$, $J$ uniform, $\varepsilon\sim\mathcal N(0,1)$).
2. $S_t=\mathcal C(\nu_\sigma(Y_t))$, $\nu_\sigma(y)=y+\sigma\xi$; $\mathcal C$ is the correction (Section 10); either may be the identity.
3. $D_t=S_t\cup(\text{reference draw})$, $\hat P_t=\Phi_h(D_t)$.

**Population (large-$N$) operator, no correction [A]:**
$$\mathcal T_{\sigma,h}(Q)=\big[(1-\lambda)(Q*N_\sigma)+\lambda P_0\big]*N_h,\qquad \hat P_t\approx\mathcal T_{\sigma,h}(\hat P_{t-1}).$$
It depends on $(\lambda,h,\sigma)$ only; it does not depend on $m$.

**Seed-conditional population operator (revision 4) [A].** Each seed has one fixed retained set with empirical law $\hat P_{\rm ref}=\frac1m\sum_j\delta_{r_j}$. Conditional on it, under the validated draw rule (every retained point used every generation), the large-$N$ operator is
$$\mathcal T^{\rm ref}_{\sigma,h}(Q)=\big[(1-\lambda)(Q*N_\sigma)+\lambda\hat P_{\rm ref}\big]*N_h.$$
Under the with-replacement rule the expectation of a reference draw is $\hat P_{\rm ref}$, so $\mathcal T^{\rm ref}$ is the mean-field operator and the per-generation draw noise is a finite-$N$ effect (Prop 3, item 4). The two operators are **not interchangeable** and define two deterministic comparators:
- **(TA) target-averaged**, $\mathcal T$ with $P_0$: averages over reference-estimation error. It is the comparator for seed-averaged summaries and the setting of Props 1–2 and of the validated tables V1–V10.
- **(SC) seed-conditional**, $\mathcal T^{\rm ref}$ with $\hat P_{\rm ref}$: keeps the seed's own reference-estimation error, so a per-seed excess over SC isolates recursive finite-sample deviation. It is the comparator for per-seed excesses ($\Delta_{R,t}$, $\Delta A_t$).
In the lineage case the distance between them is the term $(\hat\pi-\pi_3)$ of Prop 3.

---

## 4. Propositions

**Proposition 1 (contraction, fixed point) [A], numerical checks V4–V5 [G].** $\mathcal T_{\sigma,h}$ is a $(1-\lambda)$-contraction in $W_1$ (convolution is $W_1$-nonexpansive; mixing with the common component $\lambda P_0$ gives the factor $(1-\lambda)$ by coupling). Hence $W_1(\mathcal T^tQ,Q^*)\le(1-\lambda)^tW_1(Q,Q^*)$. The unique fixed point is
$$Q^*=\sum_{j\ge0}\lambda(1-\lambda)^j\,P_0*N\big(0,\,s_j^2\big),\quad s_j^2=h^2+j(h^2+\sigma^2),$$
(for $\sigma=0$, $s_j^2=(j+1)h^2$), with $\mathrm{Var}(Q^*)=v_0+c/\lambda$, $c:=h^2+(1-\lambda)\sigma^2$. Using $W_1(P_0*N_s,P_0)\le\sqrt{2/\pi}\,s$ and Jensen,
$$W_1(Q^*,P_0)\le\sqrt{\tfrac{2c}{\pi\lambda}}.$$
*Numerical status:* the variance of $Q^*$ and the bound are checked [G] (V2, V4). The full series form of $Q^*$ was **not** compared to the grid fixed point (Section 15).

*Seed-conditional extension (revision 4) [A]; checked [G]/[MC] (V13).* $\mathcal T^{\rm ref}_{\sigma,h}$ is a $(1-\lambda)$-contraction with the same proof and fixed point $Q^*_{\rm ref}=\sum_{j\ge0}\lambda(1-\lambda)^j\,\hat P_{\rm ref}*N(0,s_j^2)$. $Q^*_{\rm ref}$ and $Q^*$ are mixtures with the same weights, so by convexity of $W_1$ and nonexpansiveness of convolution
$$W_1(Q^*_{\rm ref},Q^*)\le W_1(\hat P_{\rm ref},P_0),\qquad W_1(Q^*_{\rm ref},P_0)\le W_1(\hat P_{\rm ref},P_0)+\sqrt{2c/(\pi\lambda)}.$$
A seed therefore converges toward $Q^*_{\rm ref}$, whose distance from $P_0$ contains reference-estimation error that no recursion removes.

**Proposition 2 (variance recursion and matching) [A], checked [G].** With $u_t=\mathrm{Var}(\hat P_t)-v_0$ and means stable at $\mu_0$:
$$u_t=(1-\lambda)(u_{t-1}+\sigma^2)+h^2,\quad u_0=h^2,\quad u^*=c/\lambda,\quad u_t-u^*=(1-\lambda)^t(u_0-u^*).$$
Arms with equal $c$ share $u^*$ and the bound in Proposition 1; their variance trajectories differ by $u_t^{\rm arm}-u_t^{\rm base}=-(1-\lambda)^{t+1}\sigma^2$.

**Proposition 3 (lineage-fraction dynamics) [A]; replaces revision 1.**
Setting: plain KDE recursion (noise optional, labels inherited), uniform index draws, $n_{\rm ref}=m$. Notation:
- $\pi_3=0.03$: **true target** lineage proportion (component weight).
- $\hat\pi=m_3^{\rm ref}/m$: **realised reference proportion**, where $m_3^{\rm ref}$ is the number of retained points with label 3. Under uniform retention from an iid archive, $m_3^{\rm ref}\sim\mathrm{Bin}(m,\pi_3)$.
- $\omega_t=\tfrac1N\#\{i\in D_t:\ c_i=3\}$: lineage fraction of the training set.

Conditional on $\omega_{t-1}$, the number of label-3 synthetic points is $\mathrm{Bin}(N-m,\omega_{t-1})$ exactly, so
$$\mathbb E[\omega_t\mid\omega_{t-1}]=(1-\lambda)\omega_{t-1}+\lambda\hat\pi,\qquad \mathrm{Var}[\omega_t\mid\omega_{t-1}]=\frac{(1-\lambda)\,\omega_{t-1}(1-\omega_{t-1})}{N}.$$
Consequences:
1. **Two error sources, kept separate.**
$$\omega_t-\pi_3=\underbrace{(\hat\pi-\pi_3)}_{\text{reference-estimation error (fixed per seed)}}+\underbrace{(\omega_t-\hat\pi)}_{\text{recursive deviation}}.$$
The first has variance $\pi_3(1-\pi_3)/m$ and relative sd $\sqrt{(1-\pi_3)/(\pi_3m)}$ (checked by simulation, V8 [MC]); it does not shrink with $T$ or depend on any recursion. The second has mean $(1-\lambda)^t(\mathbb E\omega_0-\hat\pi)$ [A].
2. **Stationary recursive variance (large-$N$ approximation; exact form below).** Linearising the state-dependent conditional variance at $p\approx\hat\pi$ (or $\pi_3$):
$$\mathrm{Var}(\omega_\infty\mid\hat\pi)\approx\frac{(1-\lambda)\,p(1-p)}{N\lambda(2-\lambda)}.$$
This is the large-$N$ form, not an identity, and the Gaussian/linearised form is poor near the lower boundary.

   **Exact stationary variance (revision 4) [A], checked [MC] (V11).** Write $a=1-\lambda$, $p=\hat\pi$. At stationarity the law of total variance applied to the exact conditional update, with $\mathbb E[\omega(1-\omega)\mid\hat\pi]=p(1-p)-V$ and $V:=\mathrm{Var}(\omega_\infty\mid\hat\pi)$, gives $V=\tfrac aN[p(1-p)-V]+a^2V$, hence
$$V=\frac{a\,p(1-p)}{N(1-a^2)+a}=\frac{(1-\lambda)\,\hat\pi(1-\hat\pi)}{N\lambda(2-\lambda)+(1-\lambda)}.$$
The displayed large-$N$ expression is its limit. The conditional mean is $\mathbb E[\omega_t\mid\hat\pi]=\hat\pi+(1-\lambda)^t\big(\mathbb E[\omega_0\mid\hat\pi]-\hat\pi\big)$, the expectation being over the recursion at fixed $\hat\pi$, not over archives. At stationarity the cross term between the two error sources vanishes, so the variance about $\pi_3$ is exactly $\pi_3(1-\pi_3)/m+\mathbb E[V(\hat\pi)]$.
3. **Floor and survival [A] (validated rule $n_{\rm ref}=m$ only).** $\omega_t\ge\lambda\hat\pi$ for all $t$, so lineage survival is **guaranteed** whenever $\hat\pi>0$ (under retention). Zero is absorbing only when $\hat\pi=0$, which has probability $(1-\pi_3)^m=0.97^m$ (V8). Therefore lineage *survival* is informative only for (i) non-retention arms and (ii) the stratum $\hat\pi=0$; for retained arms with $\hat\pi>0$ the informative lineage endpoint is the mass $\omega_t$, not survival.
4. **With-replacement reference rule [A]; checked [MC] (V11).** If each generation draws $n_{\rm ref}=\lambda N$ points with replacement from the retained set, the label-3 reference count is $\mathrm{Bin}(n_{\rm ref},\hat\pi)$ per generation instead of the fixed $m_3^{\rm ref}$. The conditional mean recursion is unchanged; the conditional variance gains $\lambda\hat\pi(1-\hat\pi)/N$:
$$\mathrm{Var}[\omega_t\mid\omega_{t-1}]=\frac{(1-\lambda)\omega_{t-1}(1-\omega_{t-1})+\lambda\hat\pi(1-\hat\pi)}{N},\qquad V^{\rm wr}=\frac{\hat\pi(1-\hat\pi)}{N\lambda(2-\lambda)+(1-\lambda)}=\frac{V}{1-\lambda}.$$
   The floor of item 3 no longer holds in every generation (a draw may contain no label-3 point), but zero is absorbing only when $\hat\pi=0$. The balanced-replication/subsampling rule of revision 3 is not analysed.

*Evidence status for Prop 3:* the mean recursion and conditional variance are exact consequences of inheritance and uniform draws [A]; the label chain is checked by Monte Carlo (V11) and the full KDE recursion with labels in V16 [MC]. The Monte Carlo numbers in V9 concern the **spatial** region mass, which also contains spill-over; see the caveats there.

**Proposition 4 (limits of the correction) — withdrawn in revision 4.** The identity-limit claim of revisions 2–3 was inconsistent with the noise schedule (Section 10). No analytic limit of the correction is claimed; conjectures C1–C3 (Section 10) replace it and carry no weight in any decision.
**Proposition 5 (recursion bound with correction) [A], not validated [U].** See Section 10.
**Proposition 6 (AB redundancy bound) [A]; Gaussian-pair TV inequality checked [G] (V10).** See Section 10.

---

## 5. Validation record

All entries reproduced from `validate_out.txt`. Settings are stated per table.

### V1. One-step learner bias $P_0*N_h$ (adjacent regions) [G]

Complete table (exact integration):

| $h$ | $m_{R_1}$ | $m_{R_2}$ | $m_{R_3}$ | $b_3=m_{R_3}-m_{R_3}^0$ | relative $b_3$ | $\tau(P_0*N_h)$ | $\tau/\tau_0$ |
|---|---|---|---|---|---|---|---|
| 0.3 | 0.69991 | 0.27009 | 0.03000 | −0.00000 | −0.0000 | 0.006204 | 1.369 |
| 0.6 | 0.70020 | 0.26982 | 0.02998 | −0.00002 | −0.0007 | 0.012376 | 2.730 |
| 0.9 | 0.69987 | 0.27011 | 0.03001 | +0.00001 | +0.0005 | 0.024580 | 5.422 |
| 1.2 | 0.69716 | 0.27216 | 0.03068 | +0.00068 | +0.0227 | 0.042277 | 9.326 |

Reading: one-step smoothing inflates tail mass substantially (1.4× to 9.3×). For $h\le0.9$ the region masses move by at most 0.0004 in absolute terms (the rare region by at most 0.07% relative); at $h=1.2$ the rare-region mass moves by 2.3% relative (absolute shifts up to about 0.0027, in $R_1$). Differences from the $P_0$ masses (0.69982, 0.27018, 0.03000) are arithmetic on the table [H].

### V2. Population variance recursion vs grid [G]

Grid variance vs $v_0+u_t$ over all generations ($T=600$ for $\lambda\ge0.02$, $T=1200$ for $\lambda=0.01$):

| $\lambda$ | $h$ | $\sigma$ | $\max_t\lvert\mathrm{var}_{\rm grid}-(v_0+u_t)\rvert$ | max mean drift | var(final $T$) | $v_0+c/\lambda$ |
|---|---|---|---|---|---|---|
| 0.05 | 0.60 | 0.00 | 2.90e-12 | 8.88e-16 | 15.6412 | 15.6412 |
| 0.05 | 0.60 | 0.50 | 4.84e-12 | 1.33e-15 | 20.3912 | 20.3912 |
| 0.02 | 0.30 | 0.00 | 1.78e-12 | 8.88e-16 | 12.9412 | 12.9412 |
| 0.02 | 0.30 | 0.40 | 4.80e-12 | 1.33e-15 | 20.7811 | 20.7812 |
| 0.10 | 0.90 | 0.70 | 5.07e-12 | 1.11e-15 | 20.9512 | 20.9512 |
| 0.01 | 0.30 | 0.20 | 4.85e-12 | 1.33e-15 | 21.4011 | 21.4012 |

The $10^{-4}$ gaps in two final-$T$ columns are consistent with the finite-horizon residual $(1-\lambda)^T(u_0-u^*)$ [H]; not separately tested.

### V3. Variance matching, $\lambda=0.05$, $h_0=0.6$, $T=600$ [G]

$\sigma=\kappa h_0/\sqrt{1-\lambda}$, $h=h_0\sqrt{1-\kappa^2}$.

| $\kappa$ | $\sigma$ | $h$ | var $t{=}0$ | var $t{=}5$ | var $t{=}20$ | diff vs base at $t=0,5,20$ | predicted $-(1-\lambda)^{t+1}\sigma^2$ | final var | base final var |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.1847 | 0.5724 | 8.7688 | 10.3235 | 13.1775 | −0.0324, −0.0251, −0.0116 | −0.0324, −0.0251, −0.0116 | 15.6412 | 15.6412 |
| 0.6 | 0.3694 | 0.4800 | 8.6716 | 10.2483 | 13.1427 | −0.1296, −0.1003, −0.0465 | −0.1296, −0.1003, −0.0465 | 15.6412 | 15.6412 |
| 0.9 | 0.5540 | 0.2615 | 8.5096 | 10.1229 | 13.0846 | −0.2916, −0.2256, −0.1045 | −0.2916, −0.2256, −0.1045 | 15.6412 | 15.6412 |

### V4. Population fixed points: $W_1$ vs bound; spatial endpoints [G]

$e_R=\lvert m_R-m_R^0\rvert/m_R^0$ with adjacent regions ("adj") or literal argmax ("argmax"); $e_\tau=\lvert\log(\tau/\tau_0)\rvert$; bound $=\sqrt{2c/(\pi\lambda)}$. Horizons: $T=600$ ($\lambda=0.05$), $T=1500$ ($\lambda=0.02$).

| Config / arm | $W_1^*$ | bound | $e_R^*$ (adj) | $e_R^*$ (argmax) | $e_\tau^*$ | var$^*$ |
|---|---|---|---|---|---|---|
| **$\lambda=0.05$, $h_0=0.6$** | | | | | | |
| base ($\sigma=0$) | 0.9733 | 2.1409 | 0.7358 | 0.7357 | 3.2979 | 15.6412 |
| matched $\kappa=0.3$ | 0.9719 | 2.1409 | 0.7363 | 0.7362 | 3.2961 | 15.6412 |
| matched $\kappa=0.6$ | 0.9675 | 2.1409 | 0.7379 | 0.7378 | 3.2909 | 15.6412 |
| matched $\kappa=0.9$ | 0.9596 | 2.1409 | 0.7407 | 0.7405 | 3.2822 | 15.6412 |
| unmatched $\sigma=0.5$, $h=h_0$ | 1.3648 | 2.7582 | 1.2783 | 1.2770 | 3.5737 | 20.3912 |
| **$\lambda=0.02$, $h_0=0.3$** | | | | | | |
| base ($\sigma=0$) | 0.7193 | 1.6926 | 0.4113 | 0.4112 | 3.0031 | 12.9412 |
| matched $\kappa=0.3$ | 0.7188 | 1.6926 | 0.4115 | 0.4115 | 3.0025 | 12.9412 |
| matched $\kappa=0.6$ | 0.7173 | 1.6926 | 0.4122 | 0.4122 | 3.0008 | 12.9412 |
| matched $\kappa=0.9$ | 0.7147 | 1.6926 | 0.4134 | 0.4134 | 2.9981 | 12.9412 |
| unmatched $\sigma=0.5$, $h=h_0$ | 1.7164 | 3.2655 | 1.7658 | 1.7609 | 3.7453 | 25.1912 |

Observations: the $W_1$ bound holds but is loose (about 2.2× at $\lambda=0.05$). At matched variance, population-level $W_1$ changes by at most 1.4% across $\kappa$ (0.9733→0.9596), $e_\tau$ by at most 0.5%, and $e_R$ moves by +0.7% (slightly worse); unmatched noise worsens all three. The base-arm $e_R^*\approx0.74$ and $0.41$ are **spill-over**, not rare-mode loss (Section 6).

### V5. Prop 1 contraction check, $\lambda=0.05$, $h=0.6$, $\sigma=0.3$ [G]

| $Q$ | $Q'$ | $W_1$ before | after | ratio ($\le0.95$) |
|---|---|---|---|---|
| $N(-5,1)$ | $N(6,2)$ | 11.0000 | 10.4500 | 0.9500 |
| $N(0,0.5)$ | $N(3,4)$ | 3.7600 | 3.4292 | 0.9120 |
| $N(-10,3)$ | $N(12,1)$ | 22.0000 | 20.9000 | 0.9500 |

### V6. Bandwidth-equivalence identity I1 [A] with checks [MC]

Identity [A]: $y+\sigma\xi=x_J+\sqrt{h^2+\sigma^2}\,\tilde\varepsilon$, $\tilde\varepsilon=(h\varepsilon+\sigma\xi)/\sqrt{h^2+\sigma^2}\sim\mathcal N(0,1)$.

- **(a) Coupled streams [MC; algebraic check].** $S=400$, $N=2000$, $m=100$, $T=20$, $h=0.6$, $\sigma=0.5$ ($\lambda=0.05$): $\max|D_{\rm noise}-D_{\rm bandwidth}|=1.421\times10^{-14}$.
- **(b) Independent streams [MC].** $S=1500$ seeds per arm, same settings; arm N vs arm bandwidth-coupled with independent draws:

| Statistic of $D_T$ | Noise arm mean | Bandwidth arm mean | KS $p$ |
|---|---|---|---|
| sample variance | 15.7568 | 15.8164 | 0.861 |
| count in $R_3$ ($x\ge b_{23}$) | 107.6200 | 107.8273 | 0.942 |
| sample mean | −1.3120 | −1.3178 | 0.599 |
| fraction in $R_\tau$ | 0.1319 | 0.1329 | 0.691 |

A non-significant KS test is not an equivalence test; Section 9 specifies the equivalence test.

### V7. Finite-sample vs population variance [MC]

Mean over $S=1500$ seeds of $\mathrm{var}(D_t)+h^2$ vs population $v_0+u_t$, $N=2000$, $T=20$:

| $\lambda$ | $m$ | $h$ | $\sigma$ | finite, $t=0,5,10,20$ | population, $t=0,5,10,20$ | rel. diff at 20 | SE |
|---|---|---|---|---|---|---|---|
| 0.05 | 100 | 0.6 | 0.0 | 8.79856, 10.32901, 11.49323, 13.14138 | 8.8012, 10.34854, 11.54584, 13.18916 | −0.0036 | 0.0409 |
| 0.05 | 100 | 0.6 | 0.5 | 8.79856, 11.41148, 13.40924, 16.19935 | 8.8012, 11.42308, 13.45184, 16.23635 | −0.0023 | 0.0458 |
| 0.02 | 40 | 0.3 | 0.0 | 8.52856, 8.93897, 9.31399, 9.93835 | 8.5312, 8.95491, 9.33791, 9.99705 | −0.0059 | 0.0418 |

All three differences at $t=20$ are negative (0.8–1.4 SE by arithmetic on the reported values [H]). A $(N-1)/N$ resampling-shrink origin is plausible but **not tested**.

### V8. Reference-miss probabilities [A] with simulation [MC]

$m_3^{\rm ref}\sim\mathrm{Bin}(m,\pi_3)$; "exact (A3)" uses the region mass $m_{R_3}^0=0.03000$ instead of $\pi_3$. Simulation: 40,000 draws per $m$. $\lambda=m/2000$.

| $m$ | $\lambda$ | $P(0)=0.97^m$ | exact (A3) | $P(1\text{–}2)$ exact | $P(\ge3)$ exact | sim $P(0)$ | sim $P(1\text{–}2)$ | sim $P(\ge3)$ | rel. sd formula | rel. sd sim |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.010 | 0.5438 | 0.5438 | 0.4352 | 0.0210 | 0.5407 | 0.4394 | 0.0199 | 1.271 | 1.263 |
| 40 | 0.020 | 0.2957 | 0.2957 | 0.5865 | 0.1178 | 0.2935 | 0.5878 | 0.1186 | 0.899 | 0.898 |
| 50 | 0.025 | 0.2181 | 0.2181 | 0.5927 | 0.1892 | 0.2195 | 0.5908 | 0.1898 | 0.804 | 0.803 |
| 100 | 0.050 | 0.0476 | 0.0476 | 0.3722 | 0.5802 | 0.0475 | 0.3740 | 0.5785 | 0.569 | 0.568 |
| 200 | 0.100 | 0.0023 | 0.0023 | 0.0570 | 0.9407 | 0.0022 | 0.0557 | 0.9421 | 0.402 | 0.401 |
| 400 | 0.200 | 0.0000 | 0.0000 | 0.0005 | 0.9995 | 0.0000 | 0.0004 | 0.9996 | 0.284 | 0.285 |

### V9. Spatial rare-region mass dynamics and T=20 pilot observations [MC]

**Caveat for the whole block.** These runs track the **spatial** mass of $R_3=[b_{23},\infty)$ under the KDE with $h=0.3$, not lineage. They include kernel smearing and spill-over; no lineage quantity was measured.

(a) *Stationary sd of the spatial mass $w_t=m_{R}(\hat P_t)$ vs the Prop-3 form* (within-seed sd, seeds with $m_3^{\rm ref}\ge3$, $N=2000$):

| $\lambda$ | $m$ | $T$ (burn-in) | usable seeds | empirical sd | formula | ratio |
|---|---|---|---|---|---|---|
| 0.02 | 40 | 800 (400) | **15/100** | 0.02265 | 0.01898 | 1.194 |
| 0.10 | 200 | 400 (100) | 95/100 | 0.00796 | 0.00830 | 0.959 |

Interpretation limits: the $\lambda=0.1$ row is roughly consistent with the formula. The $\lambda=0.02$ row rests on 15 usable seeds (selected by $m_3^{\rm ref}\ge3$), uses a spatial mass whose population value is inflated by spill-over at these $t$ (cf. V4: $e_R^*=0.4113$ at the fixed point for $h_0=0.3$, $\lambda=0.02$), and is **not** treated as confirmation of Prop 3.

(b) *T=20 extinction observations (legacy region-mass definition $w<\tfrac12 m_R^0=0.015$)*, $h=0.3$, $N=2000$, 1000 seeds:

| $\lambda$ | $m$ | $P(m_3^{\rm ref}=0)$ | $P(\text{extinct at }T)$ | $P(\text{ever extinct},\,t\le T)$ | mean $w_T$ | median $w_T$ |
|---|---|---|---|---|---|---|
| 0.00 | 0 | 1.000 | 0.182 | 0.317 | 0.0307 | 0.0283 |
| 0.01 | 20 | 0.531 | 0.144 | 0.257 | 0.0320 | 0.0297 |
| 0.02 | 40 | 0.287 | 0.160 | 0.274 | 0.0317 | 0.0293 |
| 0.05 | 100 | 0.044 | 0.131 | 0.257 | 0.0312 | 0.0297 |
| 0.10 | 200 | 0.004 | 0.117 | 0.223 | 0.0305 | 0.0292 |

**Treat these as pilot observations only.** They are spatial-mass extinctions at $T=20$ (not stationary for $\lambda\le0.05$, Section 7), include reference-estimation error, and do not distinguish retained from non-retained arms at the resolution of this run (0.117–0.160 with retention vs 0.182 without). The $\lambda=0$ row's $P(m_3^{\rm ref}=0)=1.000$ is trivial (no reference).

### V10. Prop 6 Gaussian-pair TV bound [G]

$\mathrm{TV}(N(0,s_0^2+\sigma^2),N(0,s_0^2))$ vs $\sigma^2/(2\sqrt2\,s_0^2)$ (numerical integration, $2.4\times10^6$ points on $[-60,60]$):

| $s_0$ | $\sigma=0.1$ | $\sigma=0.25$ | $\sigma=0.5$ | $\sigma=1.0$ |
|---|---|---|---|---|
| 0.5 | 0.0095 / 0.0141 | 0.0539 / 0.0884 | 0.1661 / 0.3536 | 0.3698 / 1.4142 |
| 1.0 | 0.0024 / 0.0035 | 0.0147 / 0.0221 | 0.0539 / 0.0884 | 0.1661 / 0.3536 |
| 2.5 | 0.0004 / 0.0006 | 0.0024 / 0.0035 | 0.0095 / 0.0141 | 0.0359 / 0.0566 |

(entries: TV / bound). Maximum TV/bound $=0.684$ ($\le1$). This checks only the Gaussian-pair inequality, **not** the claim about the correction pipeline's output laws.

### V11–V14. Revision-4 validation (scripts `validate_rev4.py`, `validate_rev4b.py`; outputs `validate_rev4_out.txt`, `validate_rev4b_out.txt`)

**V11. Lineage-chain variance, exact vs large-$N$, fixed vs with-replacement reference [MC].** Label chain only (no KDE). 40,000 independent chains per row, burn-in $\lceil\ln10^{-6}/\ln(1-\lambda)\rceil+50$ generations; Monte Carlo relative se of a variance is about 0.71% [H]. "Formula" is the exact form of Prop 3 (fixed rule) or $V^{\rm wr}$ (with-replacement); "large-$N$" is the revision-2 approximation.

| rule | $N$ | $\lambda$ | $p$ | empirical var | formula | large-$N$ | emp/formula |
|---|---|---|---|---|---|---|---|
| fixed | 2000 | 0.05 | 0.0300 | 1.4244e-04 | 1.4108e-04 | 1.4177e-04 | 1.010 |
| fixed | 2000 | 0.10 | 0.0300 | 6.8096e-05 | 6.8758e-05 | 6.8921e-05 | 0.990 |
| fixed | 100 | 0.20 | 0.0500 | 1.0408e-03 | 1.0326e-03 | 1.0556e-03 | 1.008 |
| fixed | 100 | 0.40 | 0.0500 | 4.4486e-04 | 4.4118e-04 | 4.4531e-04 | 1.008 |
| with-repl. | 2000 | 0.05 | 0.0300 | 1.4960e-04 | 1.4851e-04 | | 1.007 |
| with-repl. | 2000 | 0.10 | 0.0300 | 7.6813e-05 | 7.6398e-05 | | 1.005 |
| with-repl. | 100 | 0.20 | 0.0500 | 1.2944e-03 | 1.2908e-03 | | 1.003 |
| with-repl. | 100 | 0.40 | 0.0500 | 7.2638e-04 | 7.3529e-04 | | 0.988 |

Reading: all eight rows agree with the formula within 1.2% (about 1.7 Monte Carlo se). The simulation **cannot separate** the exact form from the large-$N$ form (they differ by at most 2.2% in these rows); the exact form rests on the derivation. The with-replacement/fixed variance ratio is $1/(1-\lambda)$ [A] (formulas 1.4851/1.4108 = 1.053; empirical 1.4960/1.4244 = 1.050 [H]). Mean transient ($N=2000$, $\lambda=0.05$, $\mathbb E\omega_0=0.015$, $\hat\pi=0.03$, $t=30$): empirical 0.026845, formula 0.026780.

**V12. Lineage-resolved population operator on the grid [G].** $dx=0.01$ on $[-30,45)$, `gaussian_filter1d`, region boundaries $b_{12},b_{23}$ of Section 2. (i) $\max\lvert\sum_kq^{(k)}_T-Q_T\rvert=1.7\times10^{-16}$ ($\lambda=0.05$, $h=0.6$, $T=59$) and $1.1\times10^{-16}$ ($\lambda=0.02$, $h=0.3$, $T=149$): the lineage-resolved operator reproduces the aggregate operator. (ii) Transition matrices $A[k,r]$ (rows: ancestry; columns: region $R_1,R_2,R_3$), target-averaged:

| $\lambda,h,T$ | $A[1,\cdot]$ | $A[2,\cdot]$ | $A[3,\cdot]$ | row sums | column sums $=m_{R_r}(Q_T)$ |
|---|---|---|---|---|---|
| 0.05, 0.6, 59 | 0.61260, 0.08388, 0.00353 | 0.04901, 0.19897, 0.02202 | 0.00018, 0.00425, 0.02557 | 0.7, 0.27, 0.03 | 0.66179, 0.28709, 0.05112 |
| 0.02, 0.3, 149 | 0.63740, 0.06163, 0.00097 | 0.03807, 0.21808, 0.01385 | 0.00005, 0.00317, 0.02678 | 0.7, 0.27, 0.03 | 0.67552, 0.28288, 0.04160 |

Row sums equal $\pi_k$ and column sums equal the region masses of $Q_T$ to the printed precision. (iii) Seed-conditional comparator, three seeds, $m=100$, $\lambda=0.05$, $h=0.6$, $T=59$: the mass of $q^{(k)}_T$ equals $\hat\pi_k+(1-\lambda)^T(\pi_k-\hat\pi_k)$ to five decimals in every case (e.g. seed 1: $\hat\pi=(0.67,0.32,0.01)$, mass $(0.67145,0.31758,0.01097)$), and $\max\lvert\sum_kq^{(k)}-Q^{\rm ref}_T\rvert\le1.4\times10^{-16}$. This validates the **population-level** identities only; no finite-sample transition matrix has been computed.

**V13. Seed-conditional fixed point and the $W_1$ bound; reference-estimation benchmark for $W_1$ [G]/[MC].** For each $m$, 40 retained-set draws; $Q^*_{\rm ref}$ from the grid operator run to the horizon of the script ($10^{-7}$ residual for $\lambda=0.05$; 800 generations for $\lambda=0.02$).

| $m$ | $\mathbb E\,W_1(\hat P_{\rm ref},P_0)$ | $\mathbb E\,W_1(Q^*_{\rm ref},Q^*)$, $\lambda=0.05$, $h=0.6$ | same, $\lambda=0.02$, $h=0.3$ |
|---|---|---|---|
| 20 | 0.7165 | 0.5356 | 0.5572 |
| 50 | 0.4515 | 0.3307 | 0.3445 |
| 100 | 0.3285 | 0.2496 | 0.2579 |
| 200 | 0.2538 | 0.1922 | 0.1991 |
| 400 | 0.1688 | 0.1342 | 0.1381 |

The inequality $W_1(Q^*_{\rm ref},Q^*)\le W_1(\hat P_{\rm ref},P_0)$ held in all 400 draws (10 cells × 40). Consequence: at $m=100$ the reference-estimation shift of the fixed point (0.25) is larger than the $W_1$ margin of Section 9.3 (0.040), so $W_1$ comparisons are interpretable only **paired** under common retained sets; The value $\mathbb E\,W_1(\hat P_{\rm ref},P_0)$ is a **benchmark** for the size of the reference-estimation discrepancy and the rationale for pairing, **not** a universal error floor: no lower bound on $W_1$ for a learner or correction is claimed, and whether a smoothed or corrected model can fall below it is not tested here. P-H2c effects are reported against C and oracle R at the same retained set.

**V14. Cross-region leakage formula vs grid [A]/[G].** The closed form $A^*[k,r]=\sum_{j}\lambda(1-\lambda)^j\pi_k[\Phi((e_{r+1}-\mu_k)/\sqrt{\sigma_k^2+s_j^2})-\Phi((e_r-\mu_k)/\sqrt{\sigma_k^2+s_j^2})]$ (finite $T$: the tail weight $(1-\lambda)^T$ at depth $T$) agrees with the grid operator to $7.4\times10^{-5}$, $6.4\times10^{-5}$, $7.4\times10^{-5}$ for $(\lambda,h,T)=(0.05,0.6,59),(0.02,0.3,149),(0.05,0.6,800)$ (grid-resolution level; the region boundary falls between grid cells). Rare-region mass by ancestry from the formula:

| $\lambda,h,T$ | lineage 1 | lineage 2 | lineage 3 | total | own-lineage retention $A[3,3]/\pi_3$ |
|---|---|---|---|---|---|
| 0.05, 0.6, 59 | 0.00354 | 0.02207 | 0.02557 | 0.05119 | 0.852 |
| 0.02, 0.3, 149 | 0.00098 | 0.01389 | 0.02678 | 0.04165 | 0.893 |
| 0.05, 0.6, 800 | 0.00424 | 0.02242 | 0.02554 | 0.05220 | 0.851 |

One-step probability of entering $R_3$ from $x=2$: $1.18\times10^{-9}$ ($h=0.6$, $\sigma=0$) and $2.25\times10^{-6}$ ($\sigma=0.5$). Cross-region mass is therefore built by diffusion over ancestry depth, not by single jumps.

### V15–V17. Revision-4 validation, pass 2 (scripts `validate_rev5_traj.py`, `validate_rev5_horizon.py`, `validate_rev5_sim.py`, `validate_rev5_sim2.py`; outputs `validate_rev5_*_out.txt`)

**V15. Endpoint trajectories and stabilisation times, uncorrected target-averaged operator, $\sigma=0$ [G].** (a) Trajectories ($dx=0.01$; $E_*$ at $t=800$ for $\lambda=0.05$, $t=1500$ for $\lambda=0.02$):

| $\lambda,h$ | $t$ | $e_R$ | $e_\tau$ | $W_1$ | $A[3,3]$ | inflow to $R_3$ |
|---|---|---|---|---|---|---|
| 0.05, 0.6 | 0 | 0.00071 | 0.99940 | 0.13618 | 0.02993 | 0.00004 |
| | 20 | 0.43554 | 3.11122 | 0.79253 | 0.02616 | 0.01690 |
| | 59 ($T_5$) | 0.70393 | 3.26934 | 0.95247 | 0.02557 | 0.02555 |
| | 118 | 0.73660 | 3.28212 | 0.97190 | 0.02554 | 0.02656 |
| | 800 | 0.73770 | 3.28256 | 0.97270 | 0.02554 | 0.02659 |
| 0.02, 0.3 | 0 | 0.00004 | 0.31200 | 0.03707 | 0.03000 | 0.00000 |
| | 20 | 0.03986 | 2.29051 | 0.42487 | 0.02874 | 0.00245 |
| | 149 ($T_5$) | 0.38674 | 2.97165 | 0.70526 | 0.02678 | 0.01482 |
| | 298 | 0.41179 | 2.98656 | 0.71826 | 0.02675 | 0.01561 |
| | 1500 | 0.41282 | 2.98709 | 0.71883 | 0.02675 | 0.01564 |

Fraction of the gap $E_0\to E_*$ closed at $T_5$: $e_R$ 0.954 ($\lambda=0.05$) and 0.937 ($\lambda=0.02$); $e_\tau$ 0.994 and 0.994; $W_1$ 0.976 and 0.980; inflow 0.961 and 0.948. (b) Sustained 5%-of-gap stabilisation time $t_{\rm stab}$ (Section 7.3 rule), $dx=0.02$, $E_*$ at residual $10^{-9}$:

| $\lambda$ | $h$ | $T_5$ | $e_R$ | $e_\tau$ | $W_1$ | $A[3,3]$ | inflow | $\max(e_R,e_\tau,W_1)/T_5$ |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 0.3 | 299 | 284 | 109 | 232 | 155 | 272 | 0.95 |
| 0.01 | 0.6 | 299 | 198 | 92 | 226 | 96 | 193 | 0.76 |
| 0.02 | 0.3 | 149 | 161 | 66 | 107 | 97 | 152 | 1.08 |
| 0.02 | 0.6 | 149 | 119 | 54 | 116 | 61 | 116 | 0.80 |
| 0.05 | 0.3 | 59 | 75 | 35 | 41 | 51 | 69 | 1.27 |
| 0.05 | 0.6 | 59 | 58 | 26 | 47 | 33 | 55 | 0.98 |
| 0.10 | 0.3 | 29 | 45 | 22 | 22 | 31 | 40 | 1.55 |
| 0.10 | 0.6 | 29 | 32 | 15 | 21 | 20 | 30 | 1.10 |
| 0.20 | 0.3 | 14 | 28 | 13 | 12 | 18 | 22 | 2.00 |
| 0.20 | 0.6 | 14 | 17 | 9 | 10 | 12 | 16 | 1.21 |

Reading: the spatial endpoint $e_R$ and the cross-lineage inflow are the slowest to settle and exceed $T_5$ in six of ten cells; the tail endpoint $e_\tau$ settles first. The two pilot configurations differ: $(0.05,0.6)$ settles at $t=58\le T_5$, $(0.02,0.3)$ at $t=161>T_5=149$. Grid note: the $e_R$ fixed-point values from different grids differ by up to about 0.003 (0.7377 at $dx=0.01$; 0.7358 in V4; 0.7400 from the closed form of V14), because the cell containing $b_{23}$ is split arbitrarily; $e_R$ read from a grid carries that uncertainty. These are properties of the **uncorrected** operator; arm B is not covered.

**V16. Full KDE lineage chain [MC].** Each chain carries positions and labels through KDE sampling ($y=x_J+h\varepsilon$, labels inherited), reference insertion and refitting, $N=2000$, $h=0.6$, one fixed retained set per row with $m_3=0.03m$ label-3 points; $\omega_t$ recorded at several late generations and the variance averaged over them; bootstrap se over chains.

| $\lambda$ | $m$ | rule | chains | var | se | var/$V_{\rm fixed}$ | var/$V_{\rm wr}$ |
|---|---|---|---|---|---|---|---|
| 0.05 | 100 | fixed | 1500 | 1.4567e-04 | 2.4% | 1.033 | 0.981 |
| 0.05 | 100 | wr | 1500 | 1.4427e-04 | 2.2% | 1.023 | 0.971 |
| 0.10 | 200 | fixed | 1500 | 6.7735e-05 | 2.0% | 0.985 | 0.887 |
| 0.10 | 200 | wr | 1500 | 7.5120e-05 | 2.0% | 1.093 | 0.983 |
| 0.20 | 400 | fixed | 1500 | 3.1728e-05 | 1.9% | 0.982 | 0.786 |
| 0.20 | 400 | wr | 1500 | 3.8426e-05 | 1.9% | 1.190 | 0.952 |
| 0.10 | 200 | fixed | 20000 | 6.8563e-05 | 0.5% | **0.997** | 0.897 |
| 0.10 | 200 | wr | 20000 | 7.7047e-05 | 0.6% | 1.121 | **1.008** |
| 0.20 | 400 | fixed | 20000 | 3.2172e-05 | 0.6% | **0.996** | 0.797 |
| 0.20 | 400 | wr | 20000 | 4.0092e-05 | 0.5% | 1.241 | **0.993** |

Reading: at $\lambda=0.05$ the two formulas differ by 5% and 1500 chains cannot separate them. At $\lambda=0.10$ and $0.20$ the fixed-rule chain matches $V_{\rm fixed}$ and the with-replacement chain matches $V_{\rm wr}$ to within 0.8%, and each rejects the other formula by many standard errors. The 4.8% shortfall of the first-pass $\lambda=0.2$ with-replacement row (2.5 se) did not persist at 20,000 chains. Because index draws ignore position, label dynamics cannot depend on the KDE positions: V16 tests the **implementation** of inheritance and reference insertion, and its agreement with V11 is expected. It does not test position-dependent endpoints.

**V17. Transition-matrix controls [MC]/[G].** $\lambda=0.05$, $h=0.6$, $T=59$, $N=2000$, $m=100$, fixed retained set (label counts 72, 27, 1, so $\hat\pi_3=0.01$); initial law $[(1-\lambda)P_0+\lambda\hat P_{\rm ref}]*N_h$; comparator $A_{\rm SC}$ from the seed-conditional grid operator. The transporter moves, in every generation, a fraction $f$ of label-1 synthetic points to $\mathcal N(8,0.6^2)$ (label kept); its prediction applies the same map inside the grid operator. 1500 chains per arm unless stated.

*Negative control (arm C, no transport):* empirical mean $A-A_{\rm SC}$ lies within $1.98$ standard errors in all nine cells (largest $\lvert z\rvert=1.98$ at $A[1,R_3]$; $\sum z^2=5.0$ over 9 cells). Row sums equal $\hat\pi_k+(1-\lambda)^T(\pi_k-\hat\pi_k)$.

*Positive control:*

| $f$ | predicted $\Delta A[1,R_3]$ | empirical (se) | share of target rare-region mass | detection $z$ | $\Delta A[3,R_3]$ (should be 0) |
|---|---|---|---|---|---|
| 0.2 | 0.54409 | 0.54456 | 1813% | 691 | 0.00000 |
| 0.002 | 0.02004 | 0.01993 (0.00019) | 67% | 103 | 0.00013 (z = 0.6) |
| 0.0005 | 0.00512 | 0.00511 (0.00013) | 17% | 40.5 | 0.00013 (z = 0.6) |

For $f=0.2$, all nine cells of the empirical transporter matrix agree with the grid prediction within $\lvert z\rvert\le1.50$ ($\sum z^2=7.8$). For $f=0.002$ and $0.0005$ the predicted row-1 changes are $(-0.01875,-0.00129,0.02004)$ and $(-0.00479,-0.00033,0.00512)$ against empirical $(-0.01837,-0.00101,0.01993)$ and $(-0.00359,-0.00053,0.00511)$; rows 2–3 are unchanged within error. Limits: these runs are unpaired with 3000 chains per arm (the $f=0.002$ and $0.0005$ rows) and test the mechanics of the measurement, not its sensitivity at $n=300$ to the effect sizes a correction will produce. Scaling the standard error by $\sqrt{10}$ for $n=300$ gives $z\approx13$ at $f=0.0005$ [H]; a $\Delta A$ below about 5% of the target rare-region mass has not been tested.

---

## 6. Rare-mode endpoints (separated)

Three endpoint families, with non-interchangeable interpretation.

### 6.1 Lineage survival and lineage mass (non-transporting arms) [U]
Defined for all arms as bookkeeping (labels inherited); informative about mode weight only for the non-transporting arms A0, C, C+N, S, DR, V. For B the label law equals that of C, so these endpoints carry no information about recovery or transport (use Sections 6.3 and 6.5).
- **Lineage mass** $\omega_{3,t}$ (weighted by the KDE weights; reference points keep their labels). Endpoints: total relative error $\epsilon^{\rm tot}_t=(\omega_{3,t}-\pi_3)/\pi_3$; decomposition $\epsilon^{\rm ref}=(\hat\pi-\pi_3)/\pi_3$ and $\epsilon^{\rm rec}_t=(\omega_{3,t}-\hat\pi)/\pi_3$ (Prop 3). Tests use absolute values.
- **Lineage survival** $1[\#\{c_i=3\}>0]$ and time to loss, reported by Kaplan–Meier. By Prop 3 this is trivially 1 under retention when $\hat\pi>0$, so it is reported (i) for A0 and (ii) within the $\hat\pi=0$ stratum, never as a win for a retained arm.

### 6.2 Conditional spatial displacement of component-3 descendants [U]
Let $G_t=\sum_{i:c_i=3}w_i\,N(x_i,h^2)\big/\sum_{i:c_i=3}w_i$ (the conditional law of label-3 descendants under the model). Endpoints:
$$D^{\rm abs}_{3,t}=W_1\big(G_t,N(\mu_3,\sigma_3^2)\big),\quad D^{\rm exc}_{3,t}=W_1\big(G_t,N(\mu_3,\sigma_3^2+h^2)\big),$$
plus mean shift $\lvert\mathbb E_{G_t}x-\mu_3\rvert$ and variance ratio $\mathrm{Var}(G_t)/\sigma_3^2$. $D^{\rm exc}$ removes the unavoidable one-step bias $h$. These are **conditional on survival**; they are paired only over seeds where both arms survive, and survivor counts are always reported. A deterministic comparator exists in principle from the lineage-resolved population operator $q^{(3)}_t=[(1-\lambda)(q^{(3)}_{t-1}*N_\sigma)+\lambda\pi_3\varphi_3]*N_h$ [A, not gridded].

### 6.3 Spatial rare-region mass with a deterministic spill-over comparator (all arms, including B) [U for comparator; G for the population trajectory itself]
- **Definition:** $m_R(Q)=Q(R_3)$, adjacent-crossing region (Section 2). Endpoint $e_R=\lvert m_R(\hat P_t)-m_R^0\rvert/m_R^0$ against $P_0$ ($m_R^0=0.03000$).
- **Deterministic spill-over comparator** $Q^{\rm pop}_t$: the grid trajectory of $\mathcal T_{\sigma,h}$ from $\hat P_0^{\rm pop}=P_0*N_h$ with the arm's $(\lambda,h,\sigma)$. Two versions (Section 3): target-averaged (TA, reference law $P_0$) for seed-averaged summaries, and seed-conditional (SC, reference law $\hat P_{\rm ref}$, $Q^{{\rm pop},{\rm ref}}_t$) for the per-seed excess $\Delta_{R,t}$. Excess $\Delta_{R,t}=m_R(\hat P_t)-m_R(Q^{\rm pop}_t)$ isolates finite-sample deviation from the no-drift prediction. The V4 values $e_R^*=0.7358$ and $0.4113$ are the comparator's own distance from $P_0$ at the fixed point.
- **For transporting correction arms (B)**, two comparators are reported and not conflated: (i) $Q^{\rm pop}_t$ of the **no-correction** operator (total effect of correction plus drift relative to pure spill-over) and (ii) a large-$N$ pseudo-population run of B with the same $(\lambda,h,\sigma,s_0,\dots)$ (B's own deterministic limit). The adequacy of (ii) requires a convergence check in $N$ [U].

### 6.4 Prohibitions and legacy
- **Spatial regional mass is never to be interpreted as component ancestry or component weight.** Evidence: at the population fixed point, $e_R^*=0.7358$ ($\lambda=0.05$, $h_0=0.6$) and $0.4113$ ($\lambda=0.02$, $h_0=0.3$) arise with no drift at all, purely from smearing of components 1–2 into $R_3$. A model that had lost component 3 entirely could therefore still exceed $\tfrac12 m_R^0$.
- **Where the spill-over comes from [A formula, G check] (V12, V14).** At $\lambda=0.05$, $h=0.6$, $T=59$ the spatial rare-region mass is $0.05119$ against $m_R^0=0.03$: lineage 3 supplies $0.02557$ (85.2% of its own mass), lineage 2 supplies $0.02207$ and lineage 1 $0.00354$. At $\lambda=0.02$, $h=0.3$, $T=149$ it is $0.04165=0.02678+0.01389+0.00098$. In the pilot the spatial error is therefore an **excess** carried by cross-lineage inflow, and an improvement in $e_R$ there means *reducing* inflow, not adding mass. The sign can reverse where the rare lineage is depleted (e.g. the $\hat\pi=0$ stratum), and $e_R$ is an absolute error, so the signed $\epsilon_R=(m_R-m_R^0)/m_R^0$ is inspected before any interpretation.
- Revision-1's $e_3$ is renamed $e_R$ and its extinction indicator $m_R<\tfrac12m_R^0$ is retired as an endpoint (kept only as the legacy observations in V9).
- Claims about mode **weight** use lineage endpoints (Section 6.1). For arm B, which can transport mass, evidence of transport is the transition matrix and conditional displacement of Section 6.5, and no mode-weight claim is made from spatial mass alone.
- Region-mass TV $\tfrac12\sum_k\lvert m_{R_k}(Q)-m_{R_k}^0\rvert$ is a secondary spatial endpoint with the same prohibition.

### 6.5 Transport evidence: ancestry-to-region transition matrix and conditional displacement (secondary, mechanistic) [U]

**Role.** This family explains *how* an arm changes the distribution. It is a secondary mechanistic endpoint: it carries no weight in the H2 decision rule, which is decided on the primary distributional endpoints (Section 12). That samples moved is not evidence that the distribution recovered.

**Transition matrix.** For the model $\hat P_t=\sum_iw_iN(x_i,h^2)$ over $D_t$ with labels $c_i$, and regions $R_1,R_2,R_3$ (Section 2):
$$A_t[k,r]=\sum_{i:\,c_i=k}w_i\,N(x_i,h^2)(R_r).$$
Identities [A]: rows sum to the lineage mass, $\sum_rA_t[k,r]=\omega_{k,t}$; columns sum to the spatial region mass, $\sum_kA_t[k,r]=m_{R_r}(\hat P_t)$. The spatial rare-region mass therefore decomposes by ancestry, $m_R=\sum_kA_t[k,3]$: own-lineage mass $A_t[3,3]$ plus inflow $I_{3,t}=\sum_{k\ne3}A_t[k,3]$. The row-normalised form $\hat A[k,r]=A[k,r]/\omega_k$ is the share of lineage $k$ sitting in region $r$.

**Deterministic comparator (spill-over removed).** The lineage-resolved population operator [A; gridded, V12]:
$$q^{(k)}_0=\pi_k\varphi_k*N_h,\qquad q^{(k)}_t=\big[(1-\lambda)(q^{(k)}_{t-1}*N_\sigma)+\lambda\pi_k\varphi_k\big]*N_h,\qquad A^{\rm pop}_t[k,r]=\int_{R_r}q^{(k)}_t.$$
By linearity $\sum_kq^{(k)}_t=Q^{\rm pop}_t$; checked on the grid to $1.7\times10^{-16}$ (V12). This is the **target-averaged** version. The **seed-conditional** version replaces $\lambda\pi_k\varphi_k$ by $\lambda\hat P^{(k)}_{\rm ref}$, the empirical reference sub-distribution carrying label $k$, whose total mass is $\hat\pi_k$ (not $\pi_k$):
$$q^{(k),{\rm ref}}_t=\big[(1-\lambda)(q^{(k),{\rm ref}}_{t-1}*N_\sigma)+\lambda\hat P^{(k)}_{\rm ref}\big]*N_h,\qquad \mathrm{mass}(q^{(k),{\rm ref}}_t)=\hat\pi_k+(1-\lambda)^t(\pi_k-\hat\pi_k).$$
The mass identity was checked on the grid for three seeds (V12). Per-seed excesses $\Delta A_t$ are taken against the seed-conditional version; seed-averaged summaries against the target-averaged version. Excess matrix $\Delta A_t=A_t-A^{\rm pop}_t$. For B two comparators are reported as in Section 6.3 (no-correction population; large-$N$ pseudo-population of B).

**Secondary summaries.** (i) own-region retention $\rho_k=A[k,k]/\omega_k$; (ii) inflow to $R_3$ and its excess $\Delta I_3=I_3-I_3^{\rm pop}$; (iii) off-diagonal excess $\lVert\Delta A\rVert_{\rm off}=\sum_{k\ne r}\lvert\Delta A[k,r]\rvert$; plus the conditional displacement $D^{\rm abs}_3$, $D^{\rm exc}_3$, mean shift and variance ratio of Section 6.2. All are reported with CIs, unadjusted, and no hypothesis test is attached.

**Why this is needed.** For B, lineage mass is blind to transport (Section 3, Prop 3): transport appears only in the **columns** of $A_t$ (position) and in displacement.

**Controls (run, V17).** *Negative:* the non-transporting arm gives $\Delta A$ consistent with zero in all nine cells (max $\lvert z\rvert=1.98$). *Positive:* an artificial transporter that moves a fraction $f$ of label-1 points into the rare-mode region reproduces the predicted $\Delta A$ in every cell (max $\lvert z\rvert=1.50$ at $f=0.2$) and is detected down to 17% of the target rare-region mass at 3000 chains per arm. The controls verify the measurement mechanics; sensitivity at $n=300$ to small effects is not established.

**Interpretation (secondary; never replaces the primary decision).** Here "H2 win" means rejection on all three primary endpoints P-H2a, P-H2b, P-H2c (Section 12.2).

| Primary endpoints (Section 12) | Transport evidence | Reading |
|---|---|---|
| H2 win | $\Delta A\approx0$ | gain from within-mode reshaping (deconvolution of the $h^2$ inflation), the regime of conjecture C3 (Section 10) |
| H2 win | $\Delta I_3\approx0$, own-lineage retention $A[3,3]$ above comparator | recovery consistent with ancestry: label-3 mass kept in or returned to its region |
| H2 win | $\Delta I_3\ne0$ in the direction that moves $m_R$ toward $m_R^0$, with no change in $A[3,3]$ | the spatial improvement comes from changed cross-lineage inflow (at the pilot an excess is reduced, Section 6.4), **not** from ancestry-preserving mode recovery; reported as such |
| no H2 win | large $\Delta A$ | movement without recovery: descriptive only, no support for H2 |
| no H2 win | $\Delta A\approx0$ | no effect |

The mechanism statement "B recovers the rare mode through mass transport" requires **both** an H2 primary win and a transport signature whose CI excludes zero; the second alone never supports it.

---

## 7. Horizon

### 7.1 Analytical residual horizon [A], arithmetic [H]
Both the $W_1$ contraction (Prop 1) and the variance transient (Prop 2) decay as $(1-\lambda)^t$. Define the 5%-residual horizon
$$T_5(\lambda)=\Big\lceil\frac{\ln0.05}{\ln(1-\lambda)}\Big\rceil.$$

| $\lambda$ | $(1-\lambda)^{20}$ [H] | $T_5(\lambda)$ [H] |
|---|---|---|
| 0.01 | 0.8179 | 299 |
| 0.02 | 0.6676 | 149 |
| 0.05 | 0.3585 | 59 |
| 0.10 | 0.1216 | 29 |
| 0.20 | 0.0115 | 14 |

[G] evidence that $T=20$ is transient: for $\lambda=0.02$, $h=0.3$, $\sigma=0$ the population variance at $t=20$ is 9.99705 vs fixed point 12.9412; for $\lambda=0.05$, $h=0.6$ it is 13.18916 vs 15.6412 (V2, V7).

### 7.2 Variance convergence does not imply distributional equilibration
The variance recursion involves only the second moment. The endpoints in Section 6 are functionals of the whole law (tail mass, region mass, lineage) with their own dynamics: the lineage mean decays at the same factor $(1-\lambda)$ [A], but its fluctuations, the tail mass and the region mass are not covered by the variance check. The validation run evaluated $e_R$ and $e_\tau$ only at fixed-point horizons (V4), not along the trajectory. Variance convergence is therefore **never** used as the stabilisation criterion for a distributional endpoint.

### 7.3 Endpoint-specific empirical stabilisation checks [U]
For each primary and secondary endpoint $E$ in each pilot configuration:
1. **Deterministic trajectory** (grid): compute $E(Q^{\rm pop}_t)$ for $t=0,\dots,\max(2T_5,T)$ and record $t_E^{(5)}$, the first $t$ such that $\lvert E_{t'}-E_*\rvert\le\max\{0.05\lvert E_*-E_0\rvert,\ \varepsilon_{\rm abs}(E)\}$ for **every** $t'\in[t,H_{\max}]$, $H_{\max}=\max(2T_5,T)$ (sustained crossing; a temporary crossing does not count; $E_*$ from the long-horizon fixed point). Edge cases, fixed in advance: if $\lvert E_*-E_0\rvert<\varepsilon_{\rm abs}(E)$ the endpoint is classed *no-drift* (stabilisation not applicable); if no such $t$ exists by $H_{\max}$, or the trajectory oscillates so that the condition fails at the end of the horizon, it is classed *not stabilised*; for finite-sample runs the same sustained rule is applied to the seed-mean trajectory. The absolute floor $\varepsilon_{\rm abs}(E)$ is endpoint-specific and not chosen here (Section 15).
2. **Finite-sample trajectory:** per-generation seed-mean of $E$ with a CI; compare successive non-overlapping windows of width 10; declare stabilised only if the window-mean change is below a tolerance $\delta_{\rm stab}$ pre-set from the pilot (fraction of the between-seed sd) **and** the late-window regression slope CI lies within $\pm\delta_{\rm stab}/10$ per generation.
3. For lineage mass compare against the analytic mean recursion; for survival use Kaplan–Meier flattening.
$\delta_{\rm stab}$ is not chosen here (Section 15).

**Separate endpoint-specific checks (decided).** Each endpoint (P-H1a/b, P-H2a/b, P-Fa/b and every secondary) gets its own stabilisation check and flag; a pass on one endpoint, or on variance, is never carried over to another. At the primary R window every endpoint is labelled *stabilised* or *not stabilised* by rules 1–2. A not-stabilised endpoint is still analysed at the R window as pre-specified, but (i) the result is labelled horizon-limited, (ii) a sensitivity analysis is run at the first window $[t_E^{\rm stab}-4,\,t_E^{\rm stab}]$ where it passes, and (iii) no equilibrium interpretation is made. Where primary and sensitivity conclusions differ, the finding is reported as horizon-dependent.

### 7.4 Two reporting horizons
Every confirmatory comparison is reported at **both**:
- **Generation-matched (G):** all arms at the same generation window (pilot: 16–20).
- **Contraction-residual-matched (R):** each configuration at the window $[T_5(\lambda)-4,\,T_5(\lambda)]$ (i.e., equal $(1-\lambda)^t\le0.05$).
Within a fixed $\lambda$ the two differ only by window location; across $\lambda$ (Section 8) they differ in meaning. The two are never pooled.

**Primary designation (decided).** Residual-matched (R) is primary for H2, AB and the factorial F; generation-matched (G) is primary for H1. The non-primary horizon is always reported as secondary. Rationale: H2, AB and F concern *accumulated* bias, whose stage is set by $(1-\lambda)^t$, so equal-residual windows match the reference-mixing exposure of the **uncorrected** operator across configurations (not a guarantee of equal convergence stage for arm B, see below); H1 is an equivalence of the same law at the same generation and does not depend on convergence stage. Cost [H]: $T_5=59$ at $\lambda=0.05$, 149 at 0.02 and 299 at 0.01, so F2 runs to 299 generations.

**Interpretation of R for arm B (revision 4).** $T_5(\lambda)=\lceil\ln0.05/\ln(1-\lambda)\rceil$ is the residual horizon of the uncorrected reference-mixing operator, and R matches that residual across configurations. Prop 5 gives arm B the bound factor $(1-\lambda)\kappa$, so B need not be at the same convergence stage as C at the same window. R is a common reference-exposure window, not an equal-convergence window. For B the estimated $\hat\kappa$ is reported, the endpoint-specific stabilisation flags of Section 7.3 apply, and when $(1-\lambda)\hat\kappa<1$ a B-specific window ending at $T_5^{B}=\lceil\ln0.05/\ln((1-\lambda)\hat\kappa)\rceil$ [H] is reported as a sensitivity analysis, never as a primary.

**Horizon evidence and proposed amendment (V15 [G]; needs ratification).** For the uncorrected target-averaged operator the sustained 5%-of-gap stabilisation time of the spatial endpoints often exceeds $T_5$: $\max(t^{\rm stab}_{e_R},t^{\rm stab}_{e_\tau},t^{\rm stab}_{W_1})/T_5$ is 0.98 at $(\lambda,h)=(0.05,0.6)$ and 1.08 at $(0.02,0.3)$, and rises to 1.27 at $(0.05,0.3)$, 1.55 at $(0.10,0.3)$ and 2.00 at $(0.20,0.3)$. $e_R$ and the cross-lineage inflow are the slowest because they accumulate by diffusion over ancestry depth (V14); $e_\tau$ settles first. $T_5$ is therefore a common exposure horizon, **not** a stabilisation horizon for the spatial endpoints. **Proposed:** $T_R(\lambda,h)=\max\{T_5(\lambda),\,t^{\rm stab}_{e_R},\,t^{\rm stab}_{e_\tau},\,t^{\rm stab}_{W_1}\}$, computed from the deterministic uncorrected trajectory and frozen before any test seed, with the R window $[T_R-4,\,T_R]$; for the two pilot configurations $T_R=59$ and $161$ (against $T_5=59$ and $149$). Until ratified, the primary R window stays at $T_5$ and the flags of Section 7.3 apply. Caveats: the criterion is relative to the gap $\lvert E_*-E_0\rvert$ and is strict when that gap is small (small $h$, large $\lambda$), where $\varepsilon_{\rm abs}$ (Section 15) would bind; the times describe the uncorrected operator, not arm B.

---

## 8. Reference-budget factorial (new)

With $N$ fixed, the proportional design $m=\lambda N$ confounds reference **information quantity** ($m$) with **mixing strength** ($\lambda$). The factorial uses the **with-replacement reference draw rule** (Section 3): every generation draws $n_{\rm ref}=\lambda N$ points with replacement from the retained set, for every cell, so $m$ and $\lambda$ vary independently at $N=2000$ under one sampling mechanism. Arm C only; CRN with nested retained sets (smaller $m$ is a prefix of larger $m$ within a seed).

| Family | What varies | What is held fixed | Pilot grid (not run) |
|---|---|---|---|
| **F1** | $m$ | $\lambda$ ($=0.05$, $n_{\rm ref}=100$) | $m\in\{20,50,100,200,400\}$ |
| **F2** | $\lambda$ | $m$ ($=100$) | $\lambda\in\{0.01,0.02,0.05,0.10\}$ ($n_{\rm ref}=20,40,100,200$) |
| **F3** | both, proportional $m=2000\lambda$ (pilot diagonal, same with-replacement rule; the validated use-every-point rule is a sensitivity analysis) | $N$ | $(m,\lambda)\in\{(20,0.01),(40,0.02),(100,0.05),(200,0.10)\}$ |

Every $m$ above appears in V8. The point $(m,\lambda)=(100,0.05)$ is shared by all three families. $\lambda=0.10$ is outside the revision-1 pilot grid and appears here (and in V8–V9) only as a comparison level; the revision-1 pilot grid is unchanged.

**What each comparison isolates**
- **F1 (vary $m$, $\lambda$ fixed):** isolates **reference information quantity**. The population operator and fixed point (Prop 1–2) are identical across $m$, as are $(1-\lambda)^t$ and $T_5$. What changes: $\mathrm{Var}(\hat\pi)=\pi_3(1-\pi_3)/m$, coverage $P(m_3^{\rm ref}=0)=0.97^m$ (V8) and reference diversity. Under the single rule the per-generation draw variance (Prop 3, item 4) depends on $(n_{\rm ref},\hat\pi)$ and is the same mechanism for every $m$; what changes is the retained set. Caveat: for small $m$ the draws contain duplicates, so the effective number of distinct reference points per generation is below $n_{\rm ref}$; information is $m$.
- **F2 (vary $\lambda$, $m$ fixed):** isolates **mixing strength**. Under CRN every $\lambda$ sees the **same** retained set (same $\hat\pi$, same coverage), so the contrast removes reference-estimation error; what changes is the contraction rate $1-\lambda$, the fixed-point offset $c/\lambda$ and $T_5(\lambda)$. Caveat: the multiplicity $n_{\rm ref}/m$ and the draw-variance term $\lambda\hat\pi(1-\hat\pi)/N$ of Prop 3 change with $\lambda$.
- **F3 (proportional):** both vary together; effects cannot be attributed. Its role is to relate the existing pilot design to F1 and F2. With $N$ fixed, the corners $\{m_{\rm lo},m_{\rm hi}\}\times\{\lambda_{\rm lo},\lambda_{\rm hi}\}$ of F1×F2 give main effects and the $m\times\lambda$ interaction; F3 is the diagonal.
- Reporting: primary at each configuration's own $T_5$ window (R, decided); fixed generation windows (G) secondary; F3 contrasts against the F1×F2 corners are descriptive.

Prop 3 predictions for F1/F2 use the with-replacement variance $V^{\rm wr}$ [A], checked on the label chain by Monte Carlo (V11); the lineage fraction through the full KDE recursion under this rule is checked (V16); other pipeline endpoints are not validated.

---

## 9. H1: noise versus bandwidth as an equivalence / negative-control experiment

**What the numerics established.** One-step smoothing inflates tails substantially without materially changing component-region masses (V1). Identity I1 [A] makes noise on synthetic samples equal in law to a generation-bandwidth change (V6). At matched variance, population-level endpoints change by at most 1.4% ($W_1$), 0.5% ($e_\tau$) and +0.7% ($e_R$) across $\kappa$ (V4); unmatched noise worsens them.

**H1 asks only:** *do the two implementations of the same law produce statistically indistinguishable outcomes?* No claim is made that Gaussian noise restores, or can restore, missing modes: an absorbed lineage (Prop 3, $\hat\pi=0$) can receive spatial mass in $R_3$ only by cross-region leakage from surviving components. The mixture has full support, so there is no empty gap, and leakage is not a single jump: a lineage-$k$ descendant at ancestry depth $j$ sits at $\mathcal N(\mu_k,\sigma_k^2+s_j^2)$, which gives the closed form of V14 [A], checked against the grid [G]. A single generation from the centre of component 2 enters $R_3$ with probability $1-\Phi\big((b_{23}-2)/\sqrt{h^2+\sigma^2}\big)=1.2\times10^{-9}$ ($h=0.6$, $\sigma=0$) [G]; the leakage in V14 is the accumulated effect over ancestry depth.

**Arms (same $\lambda$, $h$, $\sigma$, common archive/reference/$S_0$):**
- **N:** noise implementation, $y=x_J+h\varepsilon+\sigma\xi$, fit $h$.
- **H:** bandwidth-coupled implementation, $y=x_J+\sqrt{h^2+\sigma^2}\,\tilde\varepsilon$ with $\tilde\varepsilon$ drawn independently, fit $h$.
- Noise grid: variance-matched pairs $(\sigma,h(\sigma))$ at the pilot-selected $\kappa^*$ plus the stress value $\kappa=0.9$.

**Implementation check (deterministic).** Under the coupling $\tilde\varepsilon=(h\varepsilon+\sigma\xi)/\sqrt{h^2+\sigma^2}$, N and H must agree to numerical tolerance. Observed $1.421\times10^{-14}$ (V6a); pass tolerance proposed $10^{-10}$ [U]. Failure indicates a bug, not a scientific finding.

**Statistical equivalence test (paired by reference, independent recursion streams).** Within a seed, arms N and H share the archive, retained reference, $S_0$ and $D_0$ but use **independent** recursion draws (the exact coupling above is reserved for the deterministic check). The reference-estimation component of the lineage error (Prop 3) is then common to both arms and cancels in the per-seed difference. Endpoints are tested on the **signed** scale (folded absolute errors distort TOST): $\ell_\tau=\log(\tau/\tau_0)$ and $\epsilon^{\rm tot}=(\omega_3-\pi_3)/\pi_3$. Paired TOST on per-seed differences: the 90% CI of the mean difference must lie inside $(-\Delta_{{\rm eq},e},\Delta_{{\rm eq},e})$; binary endpoints use an equivalence test on discordant pairs (Newcombe-type CI inside $\pm\delta$). **Null $H_0^{(1)}$:** $\lvert\mu_N-\mu_H\rvert\ge\Delta_{{\rm eq},e}$ for at least one primary endpoint $e$; rejecting all primary nulls (intersection–union) demonstrates equivalence. Margins $\Delta_{{\rm eq},e}$: Section 9.3. Primary horizon: G.

### 9.3 Endpoint-specific equivalence margins

The universal 10% placeholder is retired. **Principle:** $\Delta_{{\rm eq},E}=\tfrac12\,\mathrm{SESOI}_E$, where $\mathrm{SESOI}_E$ is the smallest difference on endpoint $E$ that would change a conclusion in H2, AB or F. Each SESOI is anchored in a validated quantity where one exists; otherwise it is declared domain judgement for the researcher to ratify. Margins are frozen in this document; pilot data may decide **feasibility** (below) but never the margin value.

| Endpoint (signed scale) | SESOI | $\Delta_{\rm eq}$ | Anchor | Basis |
|---|---|---|---|---|
| $\ell_\tau=\log(\tau/\tau_0)$ (P-H1a) | $\ln1.25\approx0.223$ | 0.112 | V1: tail mass changes by $\ln(2.730/1.369)=0.690$ and $\ln(5.422/2.730)=0.686$ per +0.3 in one-step $h$ [H], i.e. about 0.23 per 0.1 step; SESOI is the tail effect of a 0.1 bandwidth step | anchored in V1; the "0.1 step" is judgement |
| $\epsilon^{\rm tot}=(\omega_3-\pi_3)/\pi_3$ (P-H1b) | 0.20 | 0.10 | V8: reference-estimation relative sd is at least 0.284 in every grid ($m=400$; 0.402 at $m=200$); SESOI sits below the smallest | anchored in V8; factor is judgement |
| $W_1$ (H1 secondary; primary P-H2c for H2, where no equivalence margin applies) | 0.0798 | 0.040 | $\sqrt{2/\pi}\times0.1$, the bound on the $W_1$ change from a 0.1 change in $h$ [A] | analytic |
| $\rho_v-1$ (secondary) | 0.05 | 0.025 | the pre-registered variance guard (Section 12.1) | judgement |
| $\epsilon_R=(m_R-m_R^0)/m_R^0$ (secondary) | 0.10 | 0.05 | scale only: one-step smoothing moves rare-region mass by $\le0.07\%$ for $h\le0.9$ (V1); matched allocation changes $e_R^*$ by $\le0.0049$ (V4) | judgement |
| lineage-loss proportion (secondary; A0 or $\hat\pi=0$ stratum) | 0.20 | 0.10 | none | judgement; see feasibility |

**Feasibility rule [H].** For a paired TOST with true difference 0, 90% power at $\alpha=0.05$ per side requires $\Delta_{\rm eq}\ge3.29\,\mathrm{SE}$, i.e. $\Delta_{\rm eq}\ge3.29\,s_d/\sqrt n=0.190\,s_d$ at $n=300$, where $s_d$ is the per-seed sd of the paired difference, estimated in the pilot. If an endpoint fails the rule it is classed **inconclusive by design**: it is reported with its CI and no equivalence claim is made. The margin is not widened. For binary endpoints $\Delta_{\rm eq}\ge3.29\sqrt{\psi/n}$ with $\psi$ the discordance rate; worst case $\psi=0.5$ gives $0.134$ at $n=300$, so $\Delta_{\rm eq}=0.10$ needs about 541 seeds [H] and is under-powered at 300. The $\hat\pi=0$ stratum is also small: about 14 of 300 seeds at $m=100$ ($300\times0.0476$ [H]), so lineage-loss equivalence is feasible only for A0.

### 9.4 Exploratory contrast

**Exploratory contrast (no confirmatory claim).** Variance-matched arms versus the unmatched-$h_0$ baseline compare *where smoothing is allocated* (fit vs generation) at equal variance trajectories (Prop 2), not "noise versus no noise".

---

## 10. Score estimator, correction operator and sampling dynamics (carried over; not re-validated)

Noise levels $s_1>\dots>s_L$ geometric, $s_1=s_{\max}$, $s_L=s_{\min}$. Smoothed reference law $p_s=(\sum_j\omega_j\delta_{r_j})*N_s$. Trained on $D_{\rm ref}$ only; no synthetic–original pairing.

**DSM objective [A].**
$$J(\theta)=\sum_\ell s_\ell^2\,\mathbb E_{r\sim\omega}\mathbb E_\xi\Big[\big(s_\theta(r+s_\ell\xi,s_\ell)+\xi/s_\ell\big)^2\Big].$$
**Nonparametric minimiser (S-NP) [A]:** $\hat s_s(x)=\frac1{s^2}\sum_jw_j(x)(r_j-x)$, $w_j\propto\omega_je^{-(x-r_j)^2/2s^2}$; denoiser $D_s(x)=x+s^2\hat s_s(x)$. **Neural (S-NN) [U]:** small MLP with early stopping. **Memorisation:** the unconstrained minimiser is the score of the $m$-point KDE; set $s_{\min}$ by LOO likelihood CV on the reference only.

**Correction $\mathcal C_{\ell_0}$ [U].** For input $x$: $\tilde x=x+s_0\xi$ ($s_0=s_{\ell_0}$); for $\ell=\ell_0,\dots,L$ run $K$ annealed-Langevin steps with $\alpha_\ell=\epsilon\,s_\ell^2/s_L^2$,
$$x\leftarrow x+\tfrac{\alpha_\ell}{2}s_\theta(x,s_\ell)+\sqrt{\alpha_\ell}\,z,\quad z\sim\mathcal N(0,1).$$

**Proposition 4 — withdrawn (revision 4).** *(Characterisation of C with S-NP: Section 17; C1 met, C2 partly met, C3 not met.)* Revisions 2–3 stated $s_0\le s_L\Rightarrow\mathcal C\to\mathrm{id}$. That is not available. With $s_1>\dots>s_L$ and $\ell_0\in\{1,\dots,L\}$ we have $s_0=s_{\ell_0}\ge s_L$, so the condition holds only at $\ell_0=L$. Even there the first step $\tilde x=x+s_0\xi$ perturbs the input, so $\mathcal C$ is the identity only if the Langevin updates undo the perturbation exactly, and that has not been derived. (Trivially [A]: as $\epsilon\to0$ the Langevin motion vanishes and the output is $x+s_0\xi$, not $x$.) No limit statement about $\mathcal C$ is claimed. Its displacement and mixing are to be **characterised empirically** (S-NP, as a function of $s_0,K,\epsilon$) before H2 is tested. Conjectures, all [U], used in no decision:
- **C1 (displacement):** at fixed $(K,\epsilon)$ the displacement of $\mathcal C$ increases with $s_0$.
- **C2 (reference-only regime):** for $s_0$ of the order of the mode separation and sufficient mixing, the output law depends weakly on the input, so B behaves as a reference-only sampler (to be compared with S, DR, V and oracle R).
- **C3 (within-mode regime):** for $s_0$ well below the mode separation ($d_{\min}=5$), inter-mode transport is negligible (off-diagonal $\Delta A\approx0$) and any effect is within-mode reshaping.
Each is checked with the transition matrix (Section 6.5) and displacement statistics (Section 6.2) on a pilot grid of $(s_0,K,\epsilon)$. The pilot decides which regime is operated in, not whether H2 holds.

**Proposition 5 (recursion bound with correction) [A].** If $W_1(\mathcal CQ,P_0)\le\kappa W_1(Q,P_0)+\eta$ and $(1-\lambda)\kappa<1$, with $\epsilon_x=\sqrt{2/\pi}\,x$:
$$W_1(\hat P_t,P_0)\le(1-\lambda)\big[\kappa(W_1(\hat P_{t-1},P_0)+\epsilon_\sigma)+\eta\big]+\epsilon_h,\quad d^*\le\frac{(1-\lambda)(\kappa\epsilon_\sigma+\eta)+\epsilon_h}{1-(1-\lambda)\kappa}.$$
The one-step learner bias $\epsilon_h$ is a floor for data-side correction; B can target only the accumulated excess (e.g., $h^2(1-\lambda)/\lambda$ in variance). Contraction does not require B since $(1-\lambda)<1$ already contracts. $(\kappa,\eta)$ are estimated by regression over probe models $\hat P_t$ from baseline runs (not yet run).

**Proposition 6 (AB redundancy) [A].** Per-generation TV between the correction output laws with and without prior noise $\nu_\sigma$ is at most $\sigma^2/(2\sqrt2\,s_0^2)$ (data processing plus $\mathrm{KL}\le\tfrac14(\sigma^2/s_0^2)^2$ and Pinsker); the Gaussian-pair inequality is checked in V10 [G]. If $\sigma\ll s_0$, the AB interaction is predicted near zero.

---

## 11. Baselines (equal reference budget $m$, equal $N,T$, equal compute accounting; not implemented)

- **S: stratified retention** [U]. $K=3$ k-means strata on the archive; $m_k\propto\hat\pi_k^\gamma$; weights $\omega_j=\hat\pi_{k(j)}/m_k$; uses a structural prior ($K=3$), hence an informed baseline. Lineage labels are available for retained points.
- **DR: density-ratio correction** [U]. Classifier reference vs $10^4$ synthetic draws; $\hat r=\frac{M_s}{m}\frac{c}{1-c}$ clipped at $r_{\max}$; 5-fold cross-fit; resample $N-m$ from $5(N-m)$ candidates with weights $\hat r$. Cannot create mass where $\hat P_{t-1}$ has none.
- **V: verifier selection** [U]. $v(x)=\log\Phi_{s_v}(D_{\rm ref},\omega)(x)$; threshold at the $\alpha=0.01$ LOO quantile; refill to $N-m$. Removes bias-induced outliers; cannot repair mode weights.
- **Diagnostic oracles (not used for claims):** R (reference-only KDE, no recursion), DR\* (true $p_0/\hat p$ rejection), V\* (true $p_0$ verifier).

**Arms:** A0 pure replace ($\lambda=0$); **C** uniform retention; N, H (Section 9, on top of C); S; DR; V; **B** (C + correction); **B+N** (C + noise + correction).

---

## 12. Hypotheses and endpoints (defined before the final test)

### 12.1 Endpoint definitions

Late-window average $\bar E=\tfrac15\sum_{t\in W}e(\hat P_t)$ for window $W$ = G (16–20) or R ($[T_5-4,T_5]$). **Primary window:** G for H1; R for H2, AB and F. The other window is reported as secondary (Section 7.4). Paired difference per seed $\Delta_i=\bar E_i(\text{arm})-\bar E_i(\text{comparator})$ (negative = arm better) wherever pairing under common random numbers applies.

| ID | Endpoint | Definition | Defined for |
|---|---|---|---|
| **P-H1a** | tail log-ratio (**signed** for H1) | $\ell_\tau=\log(\tau(\hat P)/\tau_0)$; $e_\tau=\lvert\ell_\tau\rvert$ elsewhere | all arms |
| **P-H1b** | lineage mass error (**signed** for H1) | $\epsilon^{\rm tot}=(\omega_3-\pi_3)/\pi_3$; absolute value elsewhere | non-transporting arms |
| **P-H2a** | rare-region mass error (spatial) | $e_R=\lvert m_R-m_R^0\rvert/m_R^0$ | all arms |
| **P-H2b** | tail log-ratio error | $e_\tau$ | all arms |
| **P-H2c** | global distributional error (revision 4) | $W_1(\hat P_t,P_0)$ | all arms |
| **P-AB** | same as P-H2a, P-H2b, P-H2c | interaction contrast | B, B+N, C, C+N |
| **P-Fa** | lineage mass error (total) | as P-H1b, with $\epsilon^{\rm ref}$, $\epsilon^{\rm rec}$ decomposition | arm C |
| **P-Fb** | tail log-ratio error | $e_\tau$ | arm C |

**Secondary mechanistic endpoints (no role in any hypothesis decision):** transition matrix $A_t$, excess $\Delta A_t$, inflow $I_3$, $\Delta I_3$, own-region retention $\rho_k$, off-diagonal excess (Section 6.5); conditional displacement (Section 6.2).

**Secondary endpoints:** mean drift; ($W_1$ is primary P-H2c for H2 and secondary elsewhere); variance ratio $\rho_v$ and $\lvert\rho_v-1\rvert$; region-mass TV; lineage survival/time-to-loss (A0 and $\hat\pi=0$ stratum only); conditional displacement $D_3^{\rm abs},D_3^{\rm exc}$; spill-over excess $\Delta_R$ against the deterministic comparator(s); stabilisation indicators; per-generation trajectories; legacy extinction indicator (V9 definition) as a descriptive item only.

**Guard (for superiority claims):** a win counts only if $\lvert\rho_v-1\rvert\le\lvert\rho_v^{C}-1\rvert+0.05$.

### 12.2 Hypotheses

**H1 (negative control / equivalence; primary horizon G).** $H_0^{(1)}$: $\lvert\mu_N-\mu_H\rvert\ge\Delta_{{\rm eq},e}$ for some $e\in\{$P-H1a, P-H1b$\}$ (non-equivalence), using the endpoint-specific margins of Section 9.3. Reject all primary nulls (IUT; paired by reference, independent recursion streams) to conclude equivalence; plus the deterministic coupling check. Own family; no alpha shared with H2.

**H2 (primary: reference-guided correction; primary horizon R).** $H_0^{(2)}$: for some baseline $b\in\{C,S,DR,V\}$, $\mathbb E[\bar E(\mathrm B)-\bar E(b)]\ge0$ on endpoint $e\in\{$P-H2a, P-H2b, P-H2c$\}$. Per endpoint, the intersection–union p-value over baselines is $p_e=\max_{b}p_{b,e}$; Holm is applied over the three endpoint-level p-values (thresholds $0.025/3$, $0.025/2$, $0.025$ in ascending order of p); per-hypothesis one-sided $\alpha=0.025$; the guard is applied. "H2 win" means rejection on all three endpoints, which under Holm is equivalent to $\max_ep_e\le0.025$ (the last Holm step), so the conjunction carries no further multiplicity penalty; partial rejections are reported with Holm-adjusted p-values.

**Decision rule.** H2 is decided on the primary distributional endpoints alone. The transition matrix and displacement (Section 6.5) are secondary mechanistic endpoints: they cannot contribute to rejection and cannot rescue a non-win. H2 must win on actual distributional recovery, not merely show that samples moved.

**Interpretation of outcomes.**
- Wins on all three: evidence that B improves the distribution relative to every specified baseline, in rare-region mass, in the tail and globally.
- Wins on P-H2a and/or P-H2b but not P-H2c: localised spatial or tail improvement; global distributional recovery is **not** established.
- Win on P-H2c only: global improvement without a rare-region or tail effect; no rare-mode claim.
- All three plus a transport signature whose CI excludes zero (Section 6.5): evidence about *how* the improvement occurred. A transport signature without the primary wins is movement without demonstrated recovery.
- Even a P-H2c win does not show that rare-component **ancestry** was recovered; that requires Sections 6.1 and 6.5. A win on spatial endpoints supports "B reduces accumulated spatial error", not by itself mode restoration (Section 6.4).

**AB (secondary interaction; gated; primary horizon R).** Tested only if H2 rejects on at least one primary endpoint. $H_0^{AB}$: $\mathbb E[\bar E(\text{B+N})]\ge\min\{\mathbb E\bar E(\text{C+N}),\mathbb E\bar E(\mathrm B)\}$ (IUT). Interaction $I=[\bar E(1,1)-\bar E(0,1)]-[\bar E(1,0)-\bar E(0,0)]$ over noise×correction at fixed variance-matched allocation; $H_0:I=0$ two-sided; $I<0$ synergy, $I\approx0$ redundancy (Prop 6 predicts $\approx0$ for $\sigma\ll s_0$). Because H1 establishes noise ≡ bandwidth, "noise" in AB means a smoothing allocation.

**F (estimation family, no confirmatory claim; primary horizon R).** Pre-specified contrasts on P-Fa, P-Fb: $C_m$ (adjacent $m$ levels at $\lambda=0.05$), $C_\lambda$ (adjacent $\lambda$ levels at $m=100$), $C_{\rm prop}$ (diagonal), with CIs; Holm within each contrast family.

---

## 13. Statistics and seeds

**Seed split.** Pilot seeds 1–100 tune $\lambda,h_0,\kappa^*,s_0,s_{\max},K,\epsilon,\gamma,r_{\max},\delta_{\rm stab},\Delta_{\rm eq}$. Test seeds 1001–1300 ($n=300$), configuration frozen, no reuse.

**Pipeline validation (must pass before hypothesis tests).**
- V-a: population variance recursion — **done** [G] (V2).
- V-b: variance matching — **done** [G] (V3).
- V-c: coupling identity — **done** [MC] (V6a).
- V-d: A0 variance $\approx v_0+(t+1)h^2$, C variance $\to v_0+h^2/\lambda$ at finite $N$, oracle ordering (R, DR\*, V\* no worse than C) — **not run** (Section 15).
- V-e: lineage inheritance and Prop 3 recursion: label chain (V11) and full KDE recursion with labels (V16) — **run**; conditional-mean transient checked on the label chain only.
- V-f: transition-matrix identities, gridded lineage-resolved comparator and negative/positive controls (Section 6.5) — **run at population level and in simulation (V12, V17)**; sensitivity at $n=300$ and realistic effect sizes not run.

**Tests.**
- Superiority endpoints: one-sided paired Wilcoxon signed-rank; paired BCa bootstrap CIs ($B=10^4$) for mean and Hodges–Lehmann differences.
- Binary endpoints: exact McNemar (paired) / Newcombe CI; equivalence of proportions for H1.
- Equivalence (H1): paired TOST on per-seed differences (pairing by shared archive, reference and $S_0$; independent recursion streams), 90% CIs, endpoint-specific margins (Section 9.3) subject to the feasibility rule.
- Transport matrix and displacement (Section 6.5): estimates with bootstrap CIs; no tests, no alpha.
- Multiplicity: Holm across the three primary endpoints of H2 (intersection–union across baselines inside each endpoint); H1 and AB as stated in Section 12.2; Holm within each contrast family of F; secondary endpoints exploratory.
- Trajectories: per-generation mean with simultaneous bootstrap bands.

**Power (normal approximation, hand calculation [H], not simulated; superiority/McNemar only, no equivalence calculation yet).** McNemar, $\delta=0.2$, discordant rate $0.4$, power 0.9: about 101 seeds ($\alpha=0.025$), 120 ($\alpha=0.0125$) and 131 ($\alpha=0.025/3$, the worst Holm level with three endpoints). Continuous: $n=300$ detects $d_z\approx0.20$ at $\alpha=0.0125$ and $d_z\approx0.21$ at $\alpha=0.025/3$ (90% power, one-sided); $n=10$ only $d_z\gtrsim1.1$. **10 seeds is inadequate** for binary endpoints; the test plan uses 300.

---

## 14. Conditions under which each hypothesis should fail

**H1 (equivalence) fails when:** the implementation of N and H differs (coupling check exceeds tolerance) — a bug; or an endpoint fails the feasibility rule of Section 9.3 (inconclusive by design, not evidence of difference); the margin is then not widened.

**H2 should fail when:**
- the retained set has no rare lineage points ($\hat\pi=0$; probability $0.97^m$): the score carries no information about component 3;
- $s_0\ll d_{\min}$ (conjecture C3): no mass transport, mode weights unchanged; a win is then a within-mode deconvolution effect;
- $s_0\gtrsim d_{\min}/2$ with ideal mixing (conjecture C2): B reduces to a reference-only sampler, matched by S, DR, V or R at equal $m$;
- $s_{\min}$ too small: memorisation (duplicates, near-zero nearest-neighbour distance to the reference);
- Langevin does not mix ($K$ or step size too small);
- $\hat\kappa$'s CI upper bound $\ge1/(1-\lambda)$: the bound in Prop 5 gives no contraction;
- reference-estimation error dominates (Prop 3): $\hat\pi$ far from $\pi_3$ caps what any retention-based method can do;
- spill-over dominates $e_R$ so that arms with equal $(\lambda,h)$ cannot be distinguished on the spatial endpoint (compare against $Q^{\rm pop}_t$);
- transport without recovery: large $\Delta A$ but no win on P-H2a/P-H2b. This is not support for H2 (Section 6.5);
- the spatial endpoint is won through a change in cross-lineage inflow ($\Delta I_3\ne0$) while own-lineage retention is unchanged: the spatial improvement is reported, but it is not evidence that the rare mode itself was recovered;
- wins on P-H2a/P-H2b without a win on P-H2c: localised improvement only (Section 12.2);
- large reference-estimation discrepancy: $\mathbb E\,W_1(\hat P_{\rm ref},P_0)=0.33$ at $m=100$ (V13) is a benchmark, not a floor; if between-seed variation in it swamps the arm contrasts, P-H2c is inconclusive unless arms are paired on common retained sets, and that pairing is mandatory.

**AB should fail (no gain over B) when:** $\sigma\ll s_0$ (Prop 6), or when either component is already at its ceiling.

**F** is uninformative at a level where $\hat\pi$ is identical across arms and the horizon is transient (generation-matched at small $\lambda$): the residual factor $(1-\lambda)^{20}$ is still 0.67–0.82 for $\lambda\le0.02$ [H].

**Global failure:** if A0 vs C is not clearly separated in the pilot at the chosen horizon, the learner is not exposing recursive bias at that horizon. Re-tune $(h_0,\lambda,T)$ or use the secondary learner (EM-GMM $K=2$; direction caveat: the KDE bias is variance inflation, opposite in sign to the shrinkage collapse in Gaussian refits).

---

## 15. Unresolved assumptions and additional numerical checks

| # | Item | Status | Check required |
|---|---|---|---|
| 1 | Lineage labels, inheritance, Prop 3 mean recursion and conditional variance | Label chain (V11) and full KDE recursion with labels (V16) checked [MC]; transient mean checked on the label chain only | Optional: transient trajectory of the mean through the full chain; per-generation survival floor under the with-replacement rule |
| 2 | Prop 3 stationary-variance approximation | Spatial-mass proxy only; $\lambda=0.02$ based on 15/100 seeds | Re-run with lineage mass, enough seeds in each $\hat\pi$ stratum |
| 3 | Lineage-resolved population operator (displacement comparator, Section 6.2) | Not gridded | Grid $q^{(3)}_t$, check mass recursion $\omega_t=(1-\lambda)\omega_{t-1}+\lambda\pi_3$ |
| 4 | Evidence of mass transport for arm B | **Decided:** ancestry-to-region transition matrix plus conditional displacement, secondary mechanistic (Section 6.5) | Items 21–22 |
| 5 | Comparator (ii) for B (large-$N$ pseudo-population) | Not run | Convergence of B's endpoints in $N$ |
| 6 | Closed-form $Q^*$ series vs grid fixed point | Variance and bound checked only | Compare $Q^*$ series to grid density ($W_1$, region masses) |
| 7 | Horizon designation | **Decided:** R for H2, AB, F; G for H1 (Section 7.4). Still open: $\delta_{\rm stab}$, window width | Pilot trajectories (Section 7.3); freeze before test seeds |
| 8 | Endpoint trajectories $e_R(t)$, $e_\tau(t)$, $W_1(t)$ on the grid | Only fixed-point values computed (V4) | Run population grid for $t\le\max(2T_5,T)$ at each pilot $(\lambda,h)$ |
| 9 | With-replacement reference rule ($n_{\rm ref}\ne m$): effect on the KDE pipeline and F1/F2 contrasts | Label-chain variance checked (V11); pipeline not implemented; revision-3 balanced-replication/subsampling rule is a sensitivity analysis only | Run the pipeline under both rules and compare with $V^{\rm wr}$ (item 29) |
| 10 | Equivalence margins | **Decided:** endpoint-specific (Section 9.3). Values proposed with anchors; those marked judgement need ratification | Pilot paired sd $s_d$ per endpoint to test feasibility ($\Delta_{\rm eq}\ge0.190\,s_d$ at $n=300$); stratum sizes for binary endpoints |
| 11 | Power figures in Section 13 | Hand normal approximation | Simulate McNemar/Wilcoxon power |
| 12 | A0 variance growth $v_0+(t+1)h^2$, finite-$N$ C convergence, oracle ordering (V-d) | Not run ($\lambda=0$ was not in the grid checks) | Run population $\lambda=0$ and finite-$N$ checks |
| 13 | Finite-$N$ variance deficit (V7, 0.8–1.4 SE, all negative) | Cause untested | Compare to $(N-1)/N$ resampling-shrink prediction with larger $S$ |
| 14 | Score estimator, correction dynamics, $(\kappa,\eta)$ estimation, Props 4–5 | Not implemented | Implement S-NP; probe-family regression; mixing diagnostics by $s_0$ |
| 15 | Baselines S, DR, V, oracles | Not implemented | Implement and unit-test (labels preserved, budget accounting) |
| 16 | Prop 6 for the full correction output | Only Gaussian-pair TV checked (V10) | Estimate TV between correction outputs with/without noise |
| 17 | Secondary learner (EM-GMM $K=2$) | Not implemented | Robustness arm |
| 18 | Stratified retention uses $K=3$ prior; $\gamma$ not chosen | Pilot | Sensitivity to strata count |
| 19 | Adaptive $\sigma_t$; 2D extension | Deferred | After 1D mechanism understood |
| 20 | Stationarity of $\lambda=0.10$ extinction at $T=20$ and of $\lambda=0.20$ cases | Pilot observations only (V9b) | Superseded by lineage endpoints at $T_5$-matched windows |
| 21 | Transition-matrix identities and gridded lineage-resolved comparator | Analytic only | Check row sums $=\omega_k$, column sums $=m_{R_r}$; check $\sum_kq^{(k)}_t$ reproduces the validated $Q^{\rm pop}_t$ trajectory |
| 22 | Sensitivity of the transition matrix to transport | Mechanics verified (V17); sensitivity at $n=300$, paired, to effects below about 5% of the target rare-region mass not run | Paired power run at $n=300$ with a weak transporter; add a within-mode transporter to test off-diagonal specificity |
| 23 | Fabrication loophole in P-H2a ($e_R$ improvable without rare-mode recovery) | **Resolved in revision 4:** $W_1$ promoted to co-primary P-H2c (Section 12.2). Residual limit: a $W_1$ win does not show ancestry recovery; and in the pilot the spatial error is an excess (V12, V14), so the loophole as first stated applied only where the rare lineage is depleted | None; see item 31 |
| 24 | Binary-endpoint feasibility for H1 (lineage loss) | $\hat\pi=0$ stratum is about 14 of 300 seeds at $m=100$; $\Delta_{\rm eq}=0.10$ needs about 541 seeds at worst-case $\psi$ [H] | Restrict to A0, enlarge $n$, or demote to descriptive |
| 25 | Pairing-by-reference design for H1 | Analytic argument only (reference-estimation error cancels) | Compare the sd of the paired difference with the unpaired sd in the pilot |
| 26 | Endpoint not stabilised at $T_5$ | Rule defined (Section 7.3), thresholds open | Choose $\delta_{\rm stab}$ in pilot; confirm sensitivity-window procedure |
| 27 | Seed-conditional vs target-averaged comparators | Operator, Prop 1 extension and population-level lineage identities checked [G]/[MC] (V12, V13); finite-sample use not run | Use SC for per-seed excess in the pilot and confirm it removes archive error from $\Delta_{R,t}$ and $\Delta A_t$ |
| 28 | Correction limits (former Prop 4) | **Withdrawn**; conjectures C1–C3 [U] | S-NP experiment over $(s_0,K,\epsilon)$ with transition matrix and displacement; derive any limit from the algorithm before restating it |
| 29 | With-replacement reference rule | Label-chain variance (V11) and full KDE chain (V16) agree with $V^{\rm wr}$ at $\lambda=0.10,0.20$; position-dependent endpoints and F1/F2 contrasts not run | Run $e_R$, $e_\tau$, $W_1$ through the KDE pipeline under both rules and compare with the SC comparator |
| 30 | Stabilisation: sustained-crossing parameters | Rule fixed in Section 7.3; $\varepsilon_{\rm abs}(E)$, $\delta_{\rm stab}$ and window width open | Choose from pilot trajectories before test seeds; confirm the no-drift and not-stabilised classes behave sensibly |
| 31 | P-H2c and the reference-estimation benchmark | $\mathbb E\,W_1(\hat P_{\rm ref},P_0)$ = 0.33 ($m=100$), 0.17 ($m=400$) [MC]; a benchmark and a pairing rationale, **not** a floor (no lower bound claimed) | Decide whether P-H2c effects are also reported as a fraction of the gap between C and oracle R; no margin is set for superiority |
| 32 | R window for arm B | Interpretation corrected (Section 7.4); $\hat\kappa$ not estimated | Estimate $\hat\kappa$ from probe runs; report the B-specific sensitivity window |
| 33 | R-horizon rule (V15) | Proposed $T_R=\max(T_5,t^{\rm stab})$ ($T_R=59$, 161 for the pilot configurations); **not applied** | **Resolved (Section 16.B.1):** rule not ratified; $T_R=T_5$ stays primary; $t^{\rm stab}$ windows are sensitivity analyses; $\varepsilon_{\rm abs}$ per Section 16.A.2 |

### Order of remaining work before freezing (from the revision-3 review; status at revision 4)

| # | Task | Status |
|---|---|---|
| 1 | $W_1$ co-primary | done in text (Section 12) |
| 2 | Conditional-reference operator | done in text; population-level checks V12, V13 |
| 3 | Correct or withdraw Prop 4 | withdrawn (Section 10) |
| 4 | Interpretation of R for B | corrected (Section 7.4) |
| 5 | Standardise F1/F2 reference sampling | done in text; lineage fraction checked through the full KDE chain (V11, V16); position-dependent endpoints not run |
| 6 | Validate lineage inheritance and Prop 3 | done for the lineage fraction (V11, V16); position-dependent endpoints not run |
| 7 | Validate transition-matrix identities and controls | population identities (V12); negative and positive controls (V17); sensitivity at $n=300$ not run |
| 8 | Grid the lineage-resolved operator | done (V12) |
| 9–12 | S-NP and correction limits; stabilisation thresholds; S, DR, V; power and equivalence feasibility | not run |


## 16. Protocol closure addendum (2026-11-09)


### A. Ratified decisions
1. **Stabilisation definition (ratified).** For endpoint E with reference level E* (the uncorrected-operator limit) and anchor E0, t_stab is the first generation t such that for every generation s ≥ t up to the horizon, |E_s − E*| ≤ max(0.05·|E* − E0|, ε_abs(E)), **and** this holds over a window of 10 consecutive generations (sustained crossing, all-subsequent condition retained).
2. **Numerical floors ε_abs (adopted as the proposed implementation).**

| Endpoint | ε_abs |
|---|---|
| W1 | 0.001 |
| Rare-region mass | 0.0001 |
| Variance ratio | 0.001 |
| Log-tail measure | 0.01 |

3. **λ = 0 arms (ratified).** T5 is undefined (no contraction; no reference mixing). **T_fix = 60** generations, generation-matched across arms; the primary terminal window is set by the experiment. T_fix is a fixed comparison horizon, not an equilibrium or convergence claim; endpoints not stabilised by generation 60 are reported as horizon-limited.
4. **Equivalence margin.** ±0.05 on late-window W1 is ratified for **Phase 2d fitting-variant comparisons only**. It is not a general W1 margin; Phase 3 margins are set separately before the test seeds.

### B. Researcher decisions on previously open items (2026-11-09)
1. **Item 33 resolved: T_R = max(T5, t_stab) is NOT ratified as the primary rule.** Primary residual-matched horizon stays **T_R = T5(λ)** (window [T5−4, T5]). Reasons: T5 gives a common reference-mixing exposure across λ; extending to the slowest endpoint gives slower configurations more recursive-sampling generations (more finite-sample accumulation); T5 of the uncorrected operator does not guarantee arm B is at the same convergence stage. Reporting structure: primary residual-matched at T5; generation-matched at the same generation across arms; stabilisation sensitivity with an endpoint-specific window ending at t_stab,E where available; arm B receives a separate B-specific stabilisation assessment (it does not replace the primary R horizon). An endpoint not stabilised at T5 stays in the primary analysis, is labelled horizon-limited, and gets the pre-specified stabilisation sensitivity analysis. No equilibrium interpretation.
2. **Phase 3 scope: the existing 1-D protocol.** 1-D three-component Gaussian mixture (rare weight 0.03), Gaussian KDE with fixed bandwidth, replace-style recursion with a fixed retained reference set and with-replacement reference draws, score-based correction (S-NP analytical; S-NN alternative), baselines S, DR, V and diagnostic oracles, H1 equivalence, H2 correction, gated AB. The 2-D accumulate-subsample / full-accumulation design is **not** part of the current confirmatory experiment; it is a deferred extension requiring its own mathematical specification and preregistration, after the 1-D correction question is answered.
3. **Seed plan retained:** pilot seeds 1–100, test seeds 1001–1300 (n = 300) (Section 1, Section 13), subject to pilot feasibility checks and design freezing before test execution. Correction to my earlier closure note: the protocol *does* contain this seed split; the only mismatch with the roadmap text was the 2-D geometry and accumulation regimes.
4. **Do not launch the 300 test seeds.** Phase 3 starts with analytical and implementation validation in this order: (1) population recursion vs numerical implementation; (2) λ = 0 control (T_fix = 60, finite-N variance); (3) with-replacement reference sampling vs lineage predictions; (4) S-NP score vs numerical differentiation of the reference KDE; (5) correction operator (Langevin updates, noise schedule, displacement, mode transport; no identity limit is assumed); (6) baselines and oracles (equal budgets, labels, compute accounting); (7) statistical feasibility (pilot variance, power, equivalence feasibility). The analytical score is not itself evidence that correction improves the distribution.

### Still open
- Item 31 (P-H2c as a fraction of the C-to-oracle-R gap) and lineage-loss equivalence feasibility.
- The R-window width (5 generations) versus the 10-generation stabilisation window are different quantities.

### D. Evidence classes for conclusions
- Phase 1 (KDE, replace-style with reinjection, 1-D): [A]/[G] for the scale-mixture identity and lineage variance; [MC] for finite-N and bandwidth results.
- Phase 2 (GMM, pre-registered): [MC], 60 paired seeds.
- Phase 2b/2c: [MC], **exploratory, post hoc**, 10–20 seeds.
- Phase 2d: [MC], pre-registered, 60 paired seeds.


Item 30 (stabilisation parameters) is resolved by A.1–A.2 (window 10, floors as tabled); Item 33 remains open.


## 17. Phase 3 validation record, steps 4–5 (2026-11-09)


Criteria: `phase3_validation_design.md` (fixed before running). Scripts `phase3_val_a.py`, `phase3_val_a2.py` (diagnostic), `phase3_val_b.py`, `phase3_val_b_analyze.py`; outputs `phase3_val_a_out.txt`, `phase3_val_a2_out.txt`, `phase3_val_b_analysis.txt`. Validation seed namespace 90001+ (pilot 1–100 and test 1001–1300 untouched). No hypothesis tested.

### Results
| Check | Outcome |
|---|---|
| V18(a) analytic score vs finite difference | **Criterion not met as written** (max rel. error 2.0e-5 vs ≤ 1e-6). Post hoc diagnostic: against a 40-digit mpmath evaluation the analytic score agrees to 4.4e-16 (s = 0.1); a 4th-order stencil gives 5e-10 for s ≥ 1 and 2e-4 only at s = 0.1. The miss is finite-difference error at small s, not the formula. [A/G] |
| V18(b) D_s = posterior mean | Met: max abs diff 5.8e-12 |
| V18(c) DSM stationarity | Met: J rises for all four perturbations (0.5730 → 0.5732–0.5740) [MC] |
| V18(d) extreme separation | Met: finite at \|x − r\|/s = 10⁵ |
| V19 s_min (LOO) | Interior optimum in 60/60 sets. Median s_min 0.79 (m = 20; range 0.25–2.46), 0.51 (m = 100), 0.35 (m = 400); all well below Silverman (1.57, 1.23, 0.92). Flag raised: s_min exceeds the planned lower s_0 in some m = 20 sets, so s_0 is taken as max(s_0, s_min) |
| V20 S-NP vs population score of P0*N_s | Descriptive. At s = 1 the mean score error in the rare region is 2.8 (m = 20), 1.3 (m = 100), 0.12 (m = 400). When the reference has no rare point (5/10, 2/10, 0/10 refs) the score at x = 8 is strongly negative (−5.3 and −5.2 at s = 1; −0.67 at s = 3, against a true −0.38 at s = 3): the score pulls samples away from the missing mode |
| V21 Langevin stationary variance | Met: sample variance within 0.71% of v/(1 − α/4v) at α/v = 0.05, 0.2, 0.5 [A/MC] |
| V22 C1 (displacement increases with s_0) | **Met**: 0 violations in 36 adjacent pairs |
| V22 C2 (output depends weakly on input, \|ΔW1\| ≤ 0.05) | **Partly met**: met for large s_0 (2.5, 6.0) at c = 0.2 (7 of 9 such cells; m = 100 and 400 at s_0 = 2.5 pass at 0.038 and 0.045), not met at s_0 ≤ 1 or at c = 0.05 (up to 0.56) |
| V22 C3 (no inter-mode transport at s_0 ≤ 1, max off-diagonal \|ΔA\| ≤ 0.005) | **Not met** in 24 of 24 cells (range 0.006–0.134). The closest cell is m = 400, s_0 = s_min, c = 0.05 (0.0059) |

### What the characterisation shows (n = 2000 inputs, 10 references per m; L = 10 levels, K = 10 steps per level)
- Fed a perfect P0 sample (input W1 to P0 = 0.08–0.09), the correction returns something at about reference quality: output W1 to P0 is 0.58–0.92 (m = 20), 0.19–0.59 (m = 100) and 0.11–0.60 (m = 400) depending on (s_0, c); at s_0 = s_min, c = 0.2 it is 0.69, 0.22, 0.13. It cannot return more than the reference holds.
- Fed the arm-B input Q* (input W1 to P0 = 0.95–0.99), it reduces W1 to 0.37–0.77 at best (s_0 ∈ {s_min, 1.0}, c = 0.2), but the output still depends on the input: at m = 400 the Q* output W1 is 0.44 against 0.13 for the same operator on a P0 input. In the pilot configuration the correction therefore recovers part of the gap, not all of it; no comparison with reference-only or with baselines has been made.
- c = 0.05 with K = 10, L = 10 is under-annealed at large s_0: at s_0 = 6 the tail mass is 0.12–0.17 (P0: 0.0045), region-3 mass 0.09–0.12 (P0: 0.030), variance 16–21 (P0: 8.4). c = 0.2 gives 0.02, 0.05–0.07 and 11–12. Choosing (L, K, c) is therefore a first-order pilot decision, not a detail.
- Inter-mode transport is present at every setting tested, increasing with s_0 and c, and decreasing with m (e.g. summed off-diagonal |ΔA| at s_0 = s_min, c = 0.05: 0.17 at m = 20, 0.04 at m = 100, 0.011 at m = 400). The correction should not be assumed to be within-mode reshaping at any s_0 ≥ s_min.

### Limits
Single-step operator on fixed input laws; the full recursion with correction was not run. One schedule family (L = 10, K = 10, geometric levels, c ∈ {0.05, 0.2}); s_0 grid relative to s_min. 10 references per m; 2000 inputs; off-diagonal ΔA has Monte Carlo noise of order 0.003, which is relevant to the 0.005 C3 threshold. C3's threshold was set by me and is not a ratified margin. S-NN, baselines, power and steps 1–3, 6–7 of the validation list are not done.


## 17.A Researcher decisions on the correction operator (2026-11-09)
| Item | Decision |
|---|---|
| C3 (no-transport regime) | Original criterion **not met** in 24 of 24 cells (threshold 0.005 was mine, unratified, with Monte Carlo noise ≈ 0.003; closest cell 0.0059). Scientific reading: the proposed no-transport regime is not supported in the tested schedule family. It is removed as an expected operating regime and is **not** replaced by a new threshold. |
| Transport in B | An expected, measurable feature, not automatically a benefit (it can move mass toward the reference's mode proportions or away from the target's). |
| Primary mechanistic diagnostic (pilot) | Full signed lineage-to-region **conditional** transition matrix ΔA = A^B − A^C (rows sum to 1, same input labels and reference set, C = no correction), with rare-region inflow/outflow, own-region retention A₃₃, and net flow J₃ = Σ_{i≠3} ω_i A_{i3} − ω₃ Σ_{j≠3} A_{3j}. The scalar off-diagonal \|ΔA\| is a summary only. |
| H2 efficacy endpoints and decision rule | Unchanged (rare-region mass, tails, W₁; transition evidence cannot reject-or-rescue H2). |
| Description of B | A reference-score-guided, annealed stochastic transport operator whose input dependence and inter-region transport depend on the entry noise scale, the Langevin step schedule and the reference sample. No mode restoration or contraction toward the target is claimed. |
| s₀ rule | s₀^actual = max(s₀^planned, s_min^LOO), computed per retained reference set and applied identically to all paired arms on that set; actual s₀ recorded, and the frequency with which the rule changes a nominal setting is reported. |
| Schedule | Selected in the pilot from a pre-specified grid (including L and K, to test whether more annealing removes the c = 0.05 tail inflation), not from the ten-reference characterisation. |
| V18(a) | Remains **"criterion not met as written"**; the 40-digit and 4th-order-stencil diagnostics are recorded separately as post hoc. |
| Independent unit | The paired seed / reference-set experiment. Inputs within a reference set estimate its conditional matrix but do not replace reference sets. Cluster/hierarchical bootstrap if both are resampled; report reference-set and input-sampling variability separately where feasible. |
| Next | Step 2, Step 6, then Step 7 (feasibility incl. H2 intersection–union, Holm, variance guard). Test seeds 1001–1300 untouched. |

## 17.B Transition-matrix read-out for V22 (arm-B input Q*, reference-set level; `phase3_val_b_transport_out.txt`)
Conditional ΔA with the no-correction comparator; one row per input label, columns R₁–R₃. Mean over 10 references per m.
- **Rare-lineage own-region retention A₃₃** (B vs C): m = 20: 0.12–0.27 vs 0.87; m = 100: 0.59–0.65 vs 0.85; m = 400: 0.75–0.84 vs 0.84 (SD across references 0.26–0.31, 0.32–0.36, 0.05–0.06). At m ≤ 100 the correction moves most rare-origin mass out of R₃ (mostly into R₂, ΔA₃₂ = +0.20 at m = 100, +0.46–0.56 at m = 20); at m = 400 it leaves most of it in place at s₀ = s_min (0.84) and moves some at s₀ = 1.0, c = 0.2 (0.75).
- **Net flow J₃** (lineage-based; Q* input): J₃^C = +0.021–0.022 at every m (non-rare lineages already sit in R₃ through smoothing leakage). J₃^B and the difference ΔJ₃ = J₃^B − J₃^C: m = 20: ≈ 0 and −0.021 to −0.024 (reference SD 0.03); m = 100: +0.011 to +0.013 and −0.009 to −0.010 (SD 0.022); m = 400: +0.022 to +0.025 and +0.001 to +0.004 (SD 0.002–0.005).
- **Reading:** at m ≤ 100 B reduces the non-rare lineage's presence in R₃ and also relocates rare-origin mass, so a change in R₃ mass cannot be read as rare-lineage preservation; at m = 400 and s₀ = s_min the matrix changes little. Reference-set variability dominates at m ≤ 100 (SD of ΔJ₃ about two to three times its mean). This is descriptive (10 references, one input law), with no inference.

## 18. Step 2 — λ = 0 control, T_fix = 60 (`phase3_step2_lam0.py`, `phase3_step2_out.txt`) [A/MC]
Fixed-horizon control, **not** an equilibrium claim. KDE (h = 0.6), N = 2000, D₀ = N fresh P₀ draws, 300 chains per stream, independent streams A (seed 91001) and B (92001).
- **Variance recursion [A]:** population variance of P̂_t is v₀ + (t+1)h² (infinite N). At finite N the sample variance of N iid draws is biased by (1 − 1/N), so E var₀(D_t) = u_t with u_t = (1 − 1/N)(u_{t−1} + h²), u₀ = (1 − 1/N)v₀, and var(P̂_t) = u_t + h².
- **Result:** stream B agrees with the finite-N prediction at every checked t within 0.23 (at most 0.8%; 0.01–0.02 for t ≤ 5). Stream A is +1.7 to +3.1 standard errors high at t = 0–5 (the t = 0 offset, +0.034, carries through t = 1–5 because the chains share their initial draws) and within 1.6 standard errors from t = 10 on. At t = 50 and 60 the population prediction is 0.53 and 0.61 above the Monte Carlo variance while the finite-N prediction is within 0.4 standard errors: the finite-N correction is real and is matched. [MC]
- **Endpoints (A/B):** W₁ 0.27/0.27 (t = 1), 0.58/0.58 (5), 0.82/0.81 (10), 1.20/1.18 (20), 1.82/1.79 (40), 2.33/2.32 (60); at t = 60 the difference between streams is 0.015 (0.7 standard errors). Region-3 mass 0.0300 → 0.108/0.104; tail mass 0.0045 → 0.256/0.259; mode count 3.0 → 4.3/4.5; no chain loses all points in R₃ over 60 generations.
- **What it shows:** with λ = 0 and a fixed-h KDE at N = 2000 the dominant failure is variance inflation by h² per generation (variance 3.4× P₀ at t = 60), with spurious tail modes, not rare-mode extinction. Caveat: mode count and tail mass at t ≥ 40 follow the inflated spread. T_fix = 60 stays a fixed horizon; nothing here says the trajectories settle.

## 19. Step 6 — arms, invariants and compute accounting (`phase3_arms.py`, `phase3_step6_out.txt`)
Arms implemented for the 1-D replace-style recursion (N = 2000, n_ref = 100 with replacement, 1900 synthetic points per generation, h = 0.6, archive 10⁴, m = 100): C, S (K = 3 k-means strata, γ = 1, weights π̂_k/m_k), DR (degree-4 logistic classifier, 5-fold cross-fit, 10⁴ synthetic, r_max = 20, 5× candidate pool), V (log-reference-KDE verifier at s_min, 1% LOO quantile threshold, refill), B (S-NP, L = 10, K = 10, c = 0.2, s₀ = max(1.0, s_min)), R (reference-only KDE), DR\* (true p₀/p̂ weights), V\* (log p₀ verifier, 1% target quantile). Seed namespace 93001+. 5 seeds, T = 5; **no endpoint comparison is reported**, to avoid selecting anything before the pilot.
- **Invariants [pass]:** |D_t| = N; reference draws all in the retained set; labels ∈ {0, 1, 2}; all values finite; retained size m = 100 with weights summing to 1, for all 8 arms.
- **Equal reference information:** the retained set is identical across C, DR, V, B, DR\*, V\*, R (same seed). S retains a different set **by design** (stratified), and is target-structure-informed (K = 3). Generation-1 synthetic input Y is identical across C, B, DR, V, DR\*, V\* (common random numbers) and differs for S because its retained set, and hence D₀, differs.
- **Compute is not equal across arms and is reported rather than matched** (per run of 5 generations, mean of 5 seeds): B 0.91 s, 9.5 × 10⁷ score evaluations, 9.5 × 10⁵ Langevin steps; DR 0.19 s, 25 classifier fits, 9.8 × 10⁴ candidate draws; DR\* 0.52 s, 9.5 × 10⁷ kernel evaluations and 4.8 × 10⁴ oracle-target evaluations; V 0.01 s, 2 × 10⁶ kernel evaluations; V\* 0.03 s; C and S negligible. Per generation B is roughly 0.18 s, so a 59-generation run is about 11 s, and 300 test seeds is ~1 h per B configuration on one core. Peak memory not yet recorded.
- **Not validated:** the DR classifier's calibration, DR\*'s weight variance, V's acceptance rate, S's weights beyond summing to 1, and any equal-compute variant. Oracles (DR\*, V\*, R) are diagnostic ceilings and use target information; they are not deployable baselines.


## 17.C Pilot schedule decision (researcher, 2026-11-09)
Use the pre-specified factorial grid s₀ ∈ {s_min, max(s_min, 1.0), max(s_min, 2.5)}, c ∈ {0.05, 0.2}, (L, K) ∈ {(10, 10), (20, 20)}, with identical or numerically indistinguishable settings removed within each reference set (at most 12 configurations per reference set). s₀ = 6 is excluded (severe tail and variance inflation in V22). (L, K) is recorded as a **joint schedule-resolution and compute factor** (100 versus 400 Langevin updates per correction), not as an isolated test of level count or of steps per level; a matched-compute diagnostic is deferred to after the pilot and only if 20 × 20 looks materially better. Schedule selection uses pre-specified numerical-validity, distributional, transport and compute criteria (Section 20), never W₁ alone. The signed conditional transition matrix is the primary mechanistic read-out, with rare-lineage retention and non-rare inflow into R₃ reported separately; uncertainty must reflect both seeds and reference sets. No efficacy claim rests on pilot selection. Test seeds 1001–1300 stay untouched until the schedule and analysis rules are locked.
Open limitations carried in Section 19: DR calibration, DR\* weight variance, V acceptance rate, equal-compute variants, peak memory.

## 20. Pilot design for schedule selection (PROPOSED — thresholds need ratification; nothing in this section has been run)

**Scope.** Pilot seeds 1–100 only. Primary configuration: λ = 0.05, m = 100 (n_ref = 100, with-replacement draws), h₀ = 0.6, N = 2000, T = T₅ = 59, R window generations 55–59; arm B vs C on the same retained set and common random numbers at generation 1 (Section 19). Other m and λ are not touched until a schedule is chosen. Mechanistic read-out per reference set: conditional ΔA, A₃₃, non-rare inflow into R₃, J₃ (Section 17.A).

**Selection procedure (in order; a configuration must pass each stage to be considered at the next).**
1. *Numerical validity:* all outputs finite in 100% of runs; implementation invariants of Section 19 pass.
2. *Distributional guardrails (proposed):* on the R window, with seed-level means: (a) variance guard as already pre-registered, |ρ_v − 1| ≤ |ρ_v^C − 1| + 0.05; (b) tail guard, tail mass τ_B ≤ 2·τ₀ (τ₀ = 0.004533), proposed as a loose sanity bound against the s₀ = 6 failure (tail mass 0.12–0.17), and not a performance target; (c) rare-region mass within [0.5, 2]·m_R⁰ (loose sanity bound). A configuration that fails a guardrail in more than 10% of seeds is dropped.
3. *Transport report:* ΔA, A₃₃, inflow, J₃ per configuration, with reference-set-level uncertainty. Reported, not thresholded (no pass/fail on transport; C3 not replaced).
4. *Compute:* updates per correction, wall-clock, peak memory (to be recorded; currently not measured).
5. *Choice among survivors:* compare P-H2a/b/c point estimates and paired uncertainty against C and R, jointly with transport and compute; prefer the cheaper schedule when 10 × 10 and 20 × 20 are not distinguishable. Ties and non-distinguishable cases are reported as such; the chosen schedule and the reason are written down before the test seeds are touched.

**Cost estimate (measured per-generation cost at m = 100, T = 59 ≈ 11 s for (10, 10)).** Six (10, 10) and six (20, 20) configurations per reference set: about 6 × 11 + 6 × 43 ≈ 320 s per seed, ≈ 8.9 core-hours for 100 seeds at one m (≈ 4.5 h on 2 cores); the baselines add little. **Proposed staging:** stage 1, seeds 1–30, all 12 configurations (~1.3 h on 2 cores), applying steps 1–2 and reporting 3–4; stage 2, seeds 31–100, only for the configurations that survive stage 1, for steps 3–5 with uncertainty. This keeps all choices inside the 100 pilot seeds. If you prefer all 12 configurations on all 100 seeds, cost is the full ≈ 4.5 h.

**Step 7 outputs from the same pilot.** Paired SDs of P-H2a/b/c (B − each of C, S, DR, V) at the chosen schedule, then power for the intersection–union comparisons, Holm over the three endpoints, and the variance guard, using the 300 test seeds; the equivalence margins for the H1 and lineage endpoints are checked for feasibility (Section 9.3).

## 20.A Pilot design decisions — FROZEN (researcher, 2026-11-09)
1. **Variance guard:** the pre-registered relative bound |ρ_v − 1| ≤ |ρ_v^C − 1| + 0.05 (a relative guardrail, not a perfect-preservation requirement).
2. **Tail guard:** KDE-smoothed tail mass ρ_tail ≤ 0.009 (a feasibility ceiling, not a target).
3. **Region-3 guard:** 0.015 ≤ ρ₃ ≤ 0.060 (sanity bound only; it is **not** a rare-mode preservation criterion, because region 3 contains non-rare ancestry; rare-lineage retention and non-rare inflow remain reported mechanistic outcomes; passing all guardrails does not establish H2 efficacy).
4. **Failure rule:** for each seed and configuration take the mean endpoint over generations 55–59 and classify the seed pass/fail separately on each guardrail. A configuration fails if more than 10% of its evaluated seeds fail **any one** guardrail (no averaging across guardrails). Individual failure rates and uncertainty intervals are reported.
5. **Stage 1:** all 12 configurations (after de-duplication) on seeds 1–30; screening and feasibility only, for obvious instability, not fine ranking; early-eliminated configurations report their stage-1 failure rate separately (not as a 100-seed estimate).
6. **Stage 2:** surviving configurations on seeds 31–100; stage-2 seeds are the primary source of schedule-comparison uncertainty; the 10% failure rule for stage-2 configurations is applied to their 100-seed results; stage-1 observations are not treated as independent confirmation.
7. **m = 100 only** for schedule selection; the selected schedule is locked and applied unchanged at all planned m in the test phase, reported separately per m, never retuned on test seeds.
8. **Signed conditional transition matrices** remain the primary mechanistic read-out; the pilot keeps three questions apart: numerical stability (finiteness, invariants, guardrails), distributional improvement (paired H2 endpoints against C and R), and ancestry change (ΔA, rare-lineage retention, non-rare inflow).
9. Paired comparisons use the same seeds, retained sets and (generation-1) input realisations across schedules and comparators; actual updates and wall-clock are recorded per configuration.
10. Step 7 (pilot paired variability → power for P-H2a/b/c under intersection–union, Holm and the variance guard) follows; 300 test seeds are not assumed sufficient.

**Operationalisation of the stage-1 screen (assistant, written before the run; open to objection).** A configuration is eliminated at stage 1 only if the lower end of the 95% Wilson interval of at least one single-guardrail failure rate exceeds 0.10 (so a narrow miss in 30 seeds does not eliminate it). Region-3 mass and tail mass are the KDE-smoothed masses of P̂_t = KDE(D_t, h); ρ_v = (var(D_t) + h²)/v₀. Transition matrices are computed on D_t with lineage labels inherited (B output inherits its input label) averaged over generations 55–59; R in the pilot is the reference-only KDE at h = 0.6 and, separately, at the LOO-selected bandwidth.

## 20.B Amendment 1 to the pilot tail guardrail (researcher, 2026-11-09; made AFTER pilot execution began)
**Audit trail.** The original frozen ceiling ρ_tail ≤ 0.009 (Section 20.A.2) was applied to the KDE-smoothed law P̂_t. After stage-1 seeds 1–2 had run, the assistant found it infeasible: the exact one-step learner P₀*N₀.₆ has tail mass 0.01238 (2.7× the unsmoothed target 0.004533), and the reference-only KDE at h = 0.6 had 0.0098–0.0102 on those two seeds, so every configuration, and the reference-only estimator, failed it. The original criterion and these two seeds' results are preserved (`pilot_stage1_audit_2seeds.pkl`; the original-threshold failure rate continues to be computed and reported alongside the amended one). This is a post-execution amendment, not an originally preregistered threshold.
**Amendment.** For this h = 0.6 pilot only, the smoothed-law tail feasibility ceiling is **ρ_tail ≤ 2 × 0.0124 = 0.0248**, calibrated to the analytic smoothing floor of P₀*N₀.₆ (fixed and reproducible, independent of the pilot's reference draws). It is not generalised to other bandwidths; the floor must be recomputed for each h.
**Conditions.** (1) The actual tail mass, the 0.0124 smoothing floor, the reference-only KDE tail mass (at h = 0.6 and at the LOO bandwidth) and the unsmoothed target 0.004533 are reported separately in all results; passing 0.0248 does not mean tail preservation was achieved (it is a feasibility guardrail). (2) Variance and region-3 guardrails, the 10% failure rule, the H2 endpoints and their comparisons are unchanged. (3) All 12 configurations continue through the planned 30 stage-1 seeds; no configuration is eliminated on the 0.009 threshold; the amended ceiling is applied consistently with the Wilson-lower-bound stage-1 screen of Section 20.A. (4) Tail mass remains a distributional check and not evidence of rare-lineage preservation; the signed conditional transition matrix remains the primary ancestry read-out.


## 21.0 FROZEN close-out of the tested S-NP Langevin schedule family (2026-11-09)
1. Stage 1 completed on pilot seeds 1–30 (m = 100) for the 12 pre-specified configurations (s₀ ∈ {s_min, max(s_min, 1.0), max(s_min, 2.5)} × c ∈ {0.05, 0.2} × (L, K) ∈ {(10, 10), (20, 20)}).
2. **0 of 12 configurations passed all feasibility guardrails** under the frozen rules (Sections 20.A and 20.B).
3. **Stage 2 was not run; test seeds 1001–1300 are untouched.**
4. The tail guard was amended (0.009 → 0.0248, Section 20.B) **after two seeds of the first, interrupted run had been observed and before any elimination decision**; the original criterion and those seeds are retained in the audit trail and the original-threshold failure rates are still reported.
5. The absolute region-3 guard was confounded by reference-set composition (reference-only KDE fails it on 33% of seeds; B's ρ₃ correlates 0.89–0.98 with it).
6. W₁ improvement over uncorrected recursion (B − C = −0.09 to −0.61) did not translate into performance matching reference-only KDE (B − R = +0.07 to +0.57).
7. **Scope of the negative result:** the tested S-NP annealed-Langevin correction schedule family and protocol at m = 100, λ = 0.05, h = 0.6, N = 2000, T = 59. It is not a rejection of reference-guided score correction in general, nor a result about other m, λ, h, score estimators, correction operators or loss functions.

## 21. Pilot stage 1 results (seeds 1–30, m = 100; screening only) — `phase3_pilot.py`, `phase3_pilot_analyze.py`, `phase3_pilot_stage1_analysis.txt`, `phase3_pilot_stage1.pkl`
Rules applied exactly as frozen in Sections 20.A and 20.B (tail ceiling 0.0248). Seeds 1–30 only; no configuration beyond stage 1 has been run; the test seeds are untouched. Evidence [MC]. A first launch was lost to a sandbox restart and rerun from scratch with identical code (resumable); the two-seed audit file `pilot_stage1_audit_2seeds.pkl` is from the rerun's first checkpoint.

**Outcome under the frozen rule: all 12 configurations are eliminated** (lower 95% Wilson bound of at least one single-guardrail failure rate exceeds 0.10 for every configuration). No stage 2 has been started.

| Guardrail | Failure-rate range over the 12 configurations | Reference-only KDE (R_h) on the same seeds | Uncorrected C |
|---|---|---|---|
| Variance (relative to C) | 0 for 11; 0.17 for s₀ = 2.5, c = 0.05, 10×10 | n/a | n/a |
| Tail ≤ 0.0248 | 0.36–1.00 (best: s₀ = 1.0, c = 0.2, 20×20 at 0.36 [0.20, 0.55]) | 0.03 | 1.00 |
| Region 3 in [0.015, 0.060] | 0.13–0.67 | 0.33 | 0.33 |
| Original tail ≤ 0.009 (audit) | 0.96–1.00 | 0.73 | 1.00 |

**What the failures are.**
- *Region-3 guard:* reference-only KDE fails it on 33% of seeds, and B's ρ₃ tracks R_h's (correlation 0.89–0.98 across seeds, 2/30 references contain no rare point). The guard therefore measures how many rare points a 100-point reference happened to hold, not the stability of the correction. A 10% failure rule cannot be met by any estimator that uses these references, including reference-only KDE.
- *Tail guard:* informative, not an artefact of the reference. R_h passes (tail 0.0129 mean, 3% failures), C fails every seed (mean 0.119), and B sits in between: mean tail 0.0205–0.0616 across configurations, about 1.75–2.0× the reference-only KDE tail at the better schedules, with the smoothing floor at 0.0124 and the unsmoothed target at 0.004533. At best B's tail is 36% failures at the 0.0248 ceiling.

**Endpoint observations (screening only; seeds 1–30; not confirmatory).**
- B reduces W₁ by 0.43–0.61 relative to C at the schedules with s₀ ≤ 1 or c = 0.2 (se ≈ 0.03–0.04); the exceptions are s₀ = 2.5, c = 0.05, 10×10 (−0.09).
- No configuration has a lower W₁ than the reference-only KDE at h = 0.6: B − R_h = +0.068 to +0.565 (R_h W₁ = 0.338, LOO-bandwidth reference-only 0.348). The best B configurations (20×20, s₀ = 1.0 or 2.5, c = 0.2) are within about +0.07 of R_h.
- Variance ratio at the R window is 1.09–1.30 for the better configurations (C: 1.80).
- Rare-lineage own-region retention A₃₃ for B is 0.11–0.73 against 0.85–0.86 for C (higher with s₀ = s_min: 0.59–0.73; 0.11–0.32 at s₀ = 1.0 or 2.5, with c = 0.2). Non-rare inflow into region 3: 0.010–0.066 for B against 0.026 for C; J₃ (lineage-based net flow): +0.002 to +0.042 for B against +0.022 for C. (Averaged over seeds where the quantity is defined; reference-set SD not shown.)
- Compute: 10×10 about 10 s per run, 20×20 about 40 s (100 and 400 updates per correction).
- s_min across the 30 references: median 0.53 (0.27–1.09). The s₀ lower-bound rule raised the nominal 1.0 in 5 seeds, and 2.5 in none; those five seeds have fewer than 3 distinct s₀ values (5 de-duplicated settings).

**Reading.** In this schedule family the correction removes most of the uncorrected recursion's error but does not beat reference-only KDE on W₁ and leaves the tail about twice as heavy as reference-only KDE. These are stage-1 observations on 30 seeds and say nothing yet about H2, which requires the confirmatory design and all baselines.


## 22. Stage-1 evidence table and existing-data mechanism analysis (no new compute; `phase3_stage1_mech.py`, `phase3_stage1_mech_out.txt`)
**Three evidence questions are kept separate: numerical stability (feasibility columns), distributional quality (W₁ columns) and ancestry mechanism (retention and inflow columns). Region-3 mass is not lineage retention.** The full 12-row table (feasibility failure rates with Wilson 95% intervals; W₁ B − C and B − R_h; own-region retention A₃₃; non-rare inflow into R₃; seconds per run; updates; seeds evaluated) is in `phase3_stage1_mech_out.txt` and not repeated here. Highlights: tail failure 36–100%, region-3 failure 13–67%, variance failure 0–17%; B − R_h = +0.068 ± 0.025 (best, s₀ = 1.0, c = 0.2, 20×20) to +0.565 ± 0.088 (worst); A₃₃ = 0.11–0.73 (C: 0.86); inflow 0.010–0.066 (C: 0.026); 10 s versus 40 s per run.

**Q1 — where does the W₁ gain come from? Not answerable at region level from stored data.** Per-generation CDFs were not stored, so W₁ cannot be split across regions without a rerun. What exists is the region mass of the sample D_t (mean over the R window), against the target (0.6998, 0.2702, 0.0300): C (0.672, 0.278, 0.050); B (smin, 0.2, 20×20) (0.678, 0.256, 0.033); B (1.0, 0.2, 20×20) (0.671, 0.254, 0.035). B largely removes C's R₃ excess (0.050 → 0.033–0.035) but R₂ ends 0.014–0.021 below the target. Reference-only region masses were not stored, so the part of this due to the 100-point reference's own composition is unidentified.

**Q2 — tail inflation: post hoc explanation, consistent with the data.** The correction ends at noise level s_min and its output is then smoothed again by the learner (h), so B's population law should resemble the reference smoothed with √(h² + s_min²), against h alone for reference-only KDE. The analytic tail mass of P₀*N_√(h²+s_min²), averaged over the 30 seeds' s_min (median 0.53, range 0.27–1.09), is 0.0237 (per-seed range 0.014–0.045), against 0.0124 for P₀*N₀.₆; the predicted variance ratio is 1.092 against 1.043. Observed at s₀ = s_min: tail 0.0239–0.0299 (best 0.0239–0.0242 for the 20×20 and c = 0.2 settings) and variance ratio 1.09–1.16. Per-seed correlation of B's tail with this prediction is +0.58 to +0.81 at s₀ = s_min. B's tail relative to R_h is 1.75–2.7× except s₀ = 2.5, c = 0.05, 10×10 (5.3×, under-annealed). Factor effects on log(tail_B / tail_R) (bootstrap over seeds): c = 0.2 vs 0.05: −0.28 [−0.31, −0.25]; 20×20 vs 10×10: −0.28 [−0.32, −0.24]; s₀ = 1.0 vs s_min: 0.00 [−0.10, +0.08]; s₀ = 2.5 vs s_min: +0.31 [+0.28, +0.35]. So the tail excess at the best schedules is mostly smoothing the schedule itself adds (about the 2× floor), not under-mixing, and the guard failures concentrate on seeds where LOO s_min is large. This is an explanation inferred after seeing the data; it has not been tested by an intervention.

**Q3 — transport.** Across the 12 configurations (Pearson, descriptive; configurations share seeds), the W₁ gap to reference-only is uncorrelated with own-region retention A₃₃ (−0.19) but tracks non-rare inflow (+0.80) and the net lineage flow J₃ (+0.99); inflow also tracks tail mass (+0.74). Within s₀ = 1.0, c = 0.2, 20×20 at the seed level, A₃₃ vs the W₁ gap is −0.33 and inflow vs the W₁ gap is +0.23 (n = 24). The lineage-3 fraction of D_t is 0.028 (C), 0.029–0.034 (B) and about 0.028 for the reference (ρ₃ proxy), i.e. lineage mass is set by the reference draws and the correction does not change it much; what changes is where lineage-labelled points sit. The W₁ gain therefore coexists with large rare-lineage relocation (A₃₃ falls to 0.11–0.73 from 0.86) and larger non-rare inflow at higher s₀. No claim that the W₁ gain is "bought" by redistribution is made: the correlations are across 12 configurations on shared seeds.

## 23. Next-question formulation (PROPOSED; nothing here has been run)
**Question:** can reference-guided score correction reduce recursive distributional degradation without introducing excessive tail distortion or cross-mode ancestry transport?
**What the stage-1 evidence says a redesign must address, not just a larger grid:** (i) the schedule's own end-of-correction smoothing doubles up with the learner's h (Q2), so tail and variance are compared to the wrong floor unless total smoothing is matched or the guards are set against the matched-smoothing floor; (ii) the region-3 guard must be conditional on the reference set; (iii) transport (non-rare inflow, J₃) co-varies with the W₁ gap to reference-only, so success has to be judged against reference-only KDE on the same reference set, not against C alone.
**Candidate mechanisms (to choose between before any run, not tested):** a correction whose final smoothing is deconvolved against the learner (matched total smoothing h_eff² = h² + s_final²); a correction with an explicit constraint on lineage-region transport; or concluding that a reference-only estimator is the target to beat for this setting.
**Design constraints:** a new section and fresh pilot seeds in a new namespace (not 1–30 and not 1001–1300); feasibility criteria and hypotheses justified and frozen before execution; no parameter sweep as such.


## 24. CDF diagnostic (seeds 1–30, previously used pilot seeds; DIAGNOSTIC, not confirmatory) — `phase3_cdf_diag.py`, `phase3_cdf_diag_analyze.py`, `phase3_cdf_diag_analysis.txt`, `phase3_cdf_diag.pkl`
Laws: C (uncorrected), B = (s₀ = max(s_min, 1.0), c = 0.2, 20×20), reference-only KDE at h (R_h), at the LOO bandwidth (R_loo) and at √(h² + s_min²) (R_match). Global CDFs, regions R₁ = (−∞, −0.105), R₂ = [−0.105, 5.58), R₃ = [5.58, ∞); regional contributions sum exactly to the total (max error 0). The rerun reproduces the stage-1 R-window W₁, tail and variance to machine precision (max difference 1.1e-16, 0, 0). The run was interrupted at 28 seeds by a sandbox restart; seeds 29–30 were completed with the same code.
- **R window, W₁ by region (R₁ | R₂ | R₃ | total):** C 0.684 | 0.235 | 0.069 | 0.988; B 0.220 | 0.149 | 0.045 | 0.414; R_h 0.160 | 0.139 | 0.039 | 0.338; R_match 0.212 | 0.144 | 0.040 | 0.396; R_loo 0.173 | 0.136 | 0.040 | 0.348.
- **Where B improves on C (paired, ±95% HW):** R₁ −0.464 ± 0.033, R₂ −0.086 ± 0.032, R₃ −0.024 ± 0.018 (total −0.574 ± 0.070): the gain is mostly in R₁ (the bulk and left tail), because C's degradation is concentrated there (signed ∫(F − F₀): C +0.27 in R₁).
- **Where B falls short of reference-only (B − R_h):** R₁ +0.060 ± 0.019, R₂ +0.010 ± 0.011, R₃ +0.006 ± 0.004, total +0.076 ± 0.024. The shortfall is concentrated in R₁, not in the rare region.
- **Against the matched-smoothing reference-only law (B − R_match):** R₁ +0.008 ± 0.007, R₂ +0.005 ± 0.009, R₃ +0.005 ± 0.004, total +0.018 ± 0.016: no detectable difference. Tail mass: B 0.0243, R_match 0.0241, R_h 0.0129 (B − R_match = +0.0001 ± 0.0004; per-seed correlation +1.00). Variance ratio: B 1.081, R_match 1.052, R_h 1.003 (+0.029 ± 0.015, correlation +0.98). Region-3 mass: 0.0312, 0.0284, 0.0283 (+0.0027 ± 0.0015). The fraction of B's excess over R_h explained by R_match is 0.99 for tail mass and 0.76 for W₁.
- **Generation-wise total W₁ (B, C, R_h):** t = 1: 0.374, 0.268, 0.338; t = 5: 0.399, 0.535, 0.338; t = 20: 0.422, 0.874, 0.338; t = 59: 0.421, 0.982, 0.338. B is worse than C at generation 1 (+0.106 ± 0.041), better from generation 5 on, and flat after about generation 5 (0.40–0.42, 0.06–0.09 above R_h throughout).
- **Reading (post hoc):** in this configuration the corrected recursion is, to within the diagnostic's resolution, a reference-only KDE with additional smoothing (the correction's terminal noise level adds to the learner's h). The W₁ improvement over C is the removal of C's own degradation; the correction adds nothing beyond a reference-only estimator at matched smoothing in tail mass and, within error, in W₁. It remains a diagnostic on seeds already used for screening.

## 25. Matched-total-smoothing hypothesis and design requirements (PROPOSED; nothing here has been run)
**Testable hypothesis H-S.** The excess tail mass (and most of the W₁ gap to R_h) of the corrected learner is attributable to the smoothing the correction adds (its terminal noise level s_L = s_min) on top of the learner's KDE smoothing h. The identity N_h ∗ N_{s_min} = N_√(h²+s_min²) establishes the equivalence only if the correction's output has that convolution structure; a Langevin stationary law is not automatically a reference KDE convolved with N_{s_min}. §24 supports H-S as a description (B ≈ R_match in tail mass, variance and W₁); it is not an identity and has not been tested by an intervention.
**Review of the three proposals received (researcher message 2026-11-10; the message ends mid-sentence at the tail guard):**
1. *Lock s₀ to s_min.* Paired comparison on the same seeds (20×20, c = 0.2): s₀ = 1.0 minus s₀ = s_min gives ΔW₁ = +0.007 ± 0.026 (n = 25) and Δtail = +0.0001 ± 0.0003; s₀ = 2.5 minus s_min gives ΔW₁ = +0.008 ± 0.022 (n = 30). The unpaired table values (+0.070 vs +0.068) used different seed sets; the paired data show no W₁ benefit from larger s₀ and a large loss of rare-lineage retention (0.59 → 0.16) and higher non-rare inflow, so the case for s₀ > s_min is not supported. Note that s₀ = s_min means **no annealing** (a constant noise level s_min over all levels), so L stops being a schedule-resolution factor and only the total Langevin evolution remains.
2. *Remove the terminal noise.* Options: (a) fit the learner at h_fit = √(max(h_min², h² − s_L²)); (b) decouple score bandwidth from thermal noise; (c) final steps at zero temperature. Cautions to settle before choosing: (c) is **not** a smoothing fix: noise-free ascent on the score of the s_min-smoothed reference moves points toward its local maxima (mode collapse and loss of within-mode spread, with possible effects on tails and on variance), and is untested; (b) changes the sampled law and needs its own stationarity analysis; (a) is the cleanest implementation of matched total smoothing but is undefined when s_L ≥ h (s_min > 0.6 occurs, up to 1.09 in 30 references) and needs a floor h_min, chosen before any run.
3. *Make both guards conditional on the same-seed reference set.* Region-3: |ρ₃(B) − ρ₃(R_h)| ≤ 0.015 is acceptable as a sanity bound (B − R_match in ρ₃ is +0.003). Tail: a bound against the matched-smoothing floor would be **vacuous for the current B** (B − R_match = +0.0001), so it cannot detect the failure H-S describes; the guard needs to be set against R_h's own tail (or the target) with a stated tolerance, and the sentence needs completing.
**Consequence for the ceiling of any redesign.** A correction that is made smoothing-matched reproduces, by construction, a reference-only KDE at the learner's bandwidth. It can restore a tie with R_h and remove the schedule-induced tail excess; it cannot exceed reference-only fidelity unless it uses information beyond the reference. Success must be defined against R_h and the best-bandwidth reference-only estimator on the same reference set, and a tie must be reported as a tie.
**Requirements for the follow-up pilot (to freeze before running):** fresh seeds in a new namespace (not 1–30 and not 1001–1300); hypotheses H-S and its competitor (score-based correction offers gain beyond matched-smoothing reference-only) stated with their tests; guards conditional on R_h; arms C, R_h, R_loo, R_match, B-original and B-matched; same retained sets and common random numbers; separate reporting of stability, distributional quality and ancestry.

## 25.A Follow-up pilot (smoothing-control) — FROZEN SPECIFICATION (written before any run; 2026-11-10)
Source: reviewer text received 2026-11-10 proposing the specification; items marked **[amended]** are the assistant's changes, made because a stated claim or criterion did not hold as written. Nothing below has been run.
**Seeds and setting.** Fresh seeds **201–230** (n = 30); seeds 1–30 are retired as diagnostic; 31–100 unused; test 1001–1300 untouched. m = 100, λ = 0.05, N = 2000, h = 0.6, T = 59, R window 55–59, identical retained sets and common random numbers across arms (generation 1 identical where the learner bandwidth is equal), s_min by LOO per reference set.
**Arms.** C (uncorrected, h = 0.6); R_h (reference-only KDE at h); R_loo (at s_min); R_match (at √(h² + s_min²)); **B_orig**: constant-level Langevin at s₀ = s_L = s_min (L = 10, K = 20 per level → 200 steps, c = 0.2), learner h = 0.6; **B_match**: identical Langevin, learner fitted at h_fit = √(max(h_min², h² − s_min²)) with **h_min = 0.15**; **R_eff**: reference-only KDE at h_eff = √(h_fit² + s_min²) = max(h, √(h_min² + s_min²)) **[amended: added as the exact comparator on clamped seeds]**. A seed is *clamped* if s_min > √(h² − h_min²) = 0.581 (8 of 30 on seeds 1–30; for clamped seeds h_eff > h and R_eff ≠ R_h; unclamped seeds have R_eff = R_h). Both strata are reported; the unclamped stratum is not used to discard seeds.
**Hypotheses.** **H-S1:** at the R window B_orig ≈ R_match: the share of B_orig's tail excess over R_h explained by R_match, (τ(R_match) − τ(R_h))/(τ(B_orig) − τ(R_h)), exceeds 0.80 (bootstrap CI over seeds), and B_orig − R_match is small in W₁. **H-S2:** B_match removes the smoothing excess: (a) conditional tail guard τ(B_match) ≤ 1.35 τ(R_h) on at least 90% of **unclamped** seeds (equivalently |τ(B_match) − τ(R_h)| ≤ 0.0045; both reported; the ratio is the primary form); (b) region-3 guard |ρ₃(B_match) − ρ₃(R_h)| ≤ 0.015 on at least 90% of unclamped seeds; (c) paired TOST equivalence of W₁(B_match) to W₁(R_eff) on all 30 seeds, margin Δ_eq = 0.040 (new, specific to this pilot; not the Phase 2d ±0.05), α = 0.05 per side; the same TOST against R_h on the unclamped stratum is secondary. **[amended: the primary comparator for W₁ on all seeds is R_eff, not R_h, because on clamped seeds B_match's exact target is R_eff.]**
**Benchmark and success bar.** R_h (and R_eff, R_loo) is the **pre-registered benchmark**, not a proven ceiling: a better reference-only estimator may exist (best bandwidth, parametric families). A tie is reported as a tie; "beats reference-only" requires W₁ superiority over R_loo and R_h with paired CIs. **[amended from "theoretical ceiling".]** Early generations (t ≤ 10) include information beyond the reference (S₀ is a fresh P₀ sample) and are transient diagnostics only, not part of the success bar.
**Ancestry.** A₃₃ and non-rare inflow are **reported, not criteria** **[amended]**: if the correction forgets its input within a generation, lineage labels carried through it no longer track mass origin, so retention is not an interpretable success measure; compared against C and B_orig on the same seeds, never against stage-1 numbers from different seeds.
**Power for the TOST (paired SD of W₁(B) − W₁(R_h) = 0.067 from §24):** n = 30: 0.90 if the true difference is 0, 0.78 at +0.01, 0.49 at +0.02; unclamped stratum n ≈ 22: 0.75, 0.64, 0.40. The stratum analyses are therefore secondary and descriptive.
**Reporting.** Stability (finiteness, guards), distributional quality (W₁ regional contributions, tail, variance, ρ₃, against R_h, R_loo, R_match, R_eff, C), ancestry (signed conditional matrices) kept separate; learner bandwidth differs between B_match (h_fit) and C/B_orig (h), so B_match − C is not a same-learner comparison and is reported as such; transient W₁ at t = 1, 5, 10, 20 recorded. Evidence tag [MC]; no claim of efficacy beyond the pre-registered contrasts.

### 25.A.R Ratification (2026-11-10, received after launch, before any 201–230 result was viewed)
User ratified all four values exactly as written in 25.A: h_min = 0.15; tail tolerance (ratio ≤ 1.35 primary, |Δτ| ≤ 0.0045 reported, unclamped seeds); TOST margin 0.040 (pilot-specific); success bar "tie with best reference-only benchmark". No amendments. Note: the run was started before ratification (assistant's error of ordering); thresholds were fixed before any result was seen. Note on margin arithmetic: 3.29×0.067/√30 = 0.0402 > 0.040, so the 0.90-power-at-zero statement is approximate (computed power 0.90 by simulation on seeds 1–30 SD; the SD is itself an estimate).

## 25.B Close-out of the smoothing-control pilot (seeds 201–230) — 2026-11-10
Analysis: `phase3_s25_analyze.py` → `phase3_s25_analysis.txt`; data `phase3_s25_pilot.pkl`; run strictly against 25.A / 25.A.R. Test seeds 1001–1300 untouched. No threshold changed after execution began.
**Audit.** (i) Spec frozen and run launched before ratification (ordering error; thresholds fixed before any result seen; ratification received before results viewed, no changes). (ii) Run interrupted at 14/30 by a sandbox process kill; resumed from the pickle (job is per-seed deterministic, SeedSequence streams), completed 30/30. (iii) Analysis bug: H-S1 bootstrap CI first computed with a generator reseeded inside the loop (degenerate interval); fixed, no criterion or threshold touched. (iv) Margin arithmetic: 3.29×0.067/√30 = 0.0402 > 0.040, so "0.90 power at zero difference" is approximate (SD is itself an estimate); immaterial because H-S2c failed by a wide margin. (v) Strata on these seeds: 7 clamped, 23 unclamped (8 and 22 on seeds 1–30).
**Results [MC].** R-window means (W1 / tau / rho3 / rv): C 1.033/0.120/0.051/1.78; R_h 0.312/0.0125/0.0297/1.014; R_loo 0.317/0.0119/0.0297/1.012; R_match 0.363/0.0215/0.0298/1.054; R_eff 0.327/0.0151/0.0297/1.026; B_orig 0.386/0.0210/0.0292/1.048; B_match 0.355/0.0148/0.0284/1.008.
- **H-S1: supported for tail.** (τ(R_match)−τ(R_h))/(τ(B_orig)−τ(R_h)) = 1.06 (bootstrap CI 1.02–1.12 > 0.80); per-seed tail correlation 0.995. The W1 clause ("small") is only partly met: B_orig − R_match = +0.023 (95% CI 0.004–0.042).
- **H-S2a: pass.** Unclamped: mean τ ratio 0.98, 100% ≤ 1.35 (and 100% within 0.0045). B_orig on the same seeds: mean ratio 1.52, 43% ≤ 1.35. (The 1.62 figure is from seeds 1–30.) Clamped seeds (outside the guard's scope; reported): mean ratio 2.22, 14% pass.
- **H-S2b: fail.** |ρ3(B_match) − ρ3(R_h)| ≤ 0.015 on 78% of unclamped seeds (< 90%); mean |Δ| 0.0099; mean ρ3 is lower than R_h's (0.0284 vs 0.0297).
- **H-S2c: fail.** W1(B_match) − W1(R_eff) = +0.027, 90% CI [0.006, 0.049]; upper bound exceeds 0.040, so equivalence is not shown; the 95% CI [0.003, 0.053] excludes 0, so B_match is also worse. Secondary, unclamped vs R_h: +0.038, 90% CI [0.011, 0.065], not equivalent.
- **Superiority/tie bar: not met.** B_match − R_h = +0.042 (CI 0.016–0.068); − R_loo = +0.038 (0.013–0.064); − R_eff = +0.027 (0.003–0.053). No claim of efficacy against reference-only.
- **vs C:** B_orig − C = −0.647 (same learner); B_match − C = −0.678 (different learner; not a same-learner comparison).
- **Transient W1 (t = 1, 5, 10, 20):** C 0.260/0.523/0.684/0.851; B_orig 0.294/0.351/0.378/0.405; B_match 0.246/0.309/0.343/0.357; reference-only R_h is constant 0.312 per seed. B_match is below R_h at t=1 (S_0 still holds 1900 fresh P_0 draws), about equal at t=5, above from t=10, and plateaus near the R-window value 0.355 by t=20 (it does not keep rising). A_33: C 0.85, B_orig 0.73, B_match 0.78 (reported only).
**Conclusion.** The matched-smoothing hypothesis is confirmed by intervention for tail mass and variance on unclamped seeds, and the post-hoc §22 account (B ≈ R_match) is now intervention-supported for those two quantities [MC]. The stronger claim, that B_match ties its exact reference-only target in W1 and in region-3 mass, is refuted at this design's power. Net for the S-NP Langevin line in the 1-D stationary setting (m = 100, λ = 0.05): partial/negative finding — it removes most uncorrected degradation but does not match or beat the reference-only benchmark; clamped seeds (7/30) still fail the tail guard. Claim scope: the tested constant-level s_min Langevin (200 steps, c = 0.2) with learner h_fit; not score correction in general.
**Residual W1 gap, not explained [U/H].** Candidate explanations are hypotheses only and were not tested: (1) incomplete cross-region mixing within a generation shifting regional masses (the ρ3 shortfall is consistent in sign but untested; no mixing diagnostic was run); (2) ULA discretisation bias (V21 measured the stationary-variance bias at 0.71%, which on its face is too small to produce +0.027 W1); (3) h_fit < h changing the learner's own resampling smoothness. The earlier statement that 400 steps erase S_0 at t = 1 is not supported: at s_min with 200 steps B_match already beats C at t = 1.
**Status.** Closed. No further S-NP Langevin parameter sweep. Any further question (e.g. a mixing diagnostic for the residual) needs fresh seeds (new namespace) and pre-frozen criteria. Evidence tag [MC]; distinct questions (stability, distributional quality, ancestry) kept separate.
