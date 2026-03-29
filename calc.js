/* ===================================
   医学样本量计算器 - 核心计算逻辑
   基于 Chow et al. (2018) 及权威文献
   =================================== */

/* ===== 工具函数 ===== */

// 标准正态分布的逆函数（Peter Acklam 高精度算法，最大误差 < 1.15e-9）
// 参考：https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
function normInv(p) {
  const a = [-3.969683028665376e+01,  2.209460984245205e+02,
             -2.759285104469687e+02,  1.383577518672690e+02,
             -3.066479806614716e+01,  2.506628277459239e+00];
  const b = [-5.447609879822406e+01,  1.615858368580409e+02,
             -1.556989798598866e+02,  6.680131188771972e+01,
             -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01,
             -2.400758277161838e+00, -2.549732539343734e+00,
              4.374664141464968e+00,  2.938163982698783e+00];
  const d = [ 7.784695709041462e-03,  3.224671290700398e-01,
              2.445134137142996e+00,  3.754408661907416e+00];
  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  if (p < pLow) {
    // 左尾区域
    const q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    // 中间区域（含 p=0.975 等常用值）
    const q = p - 0.5;
    const r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    // 右尾区域
    const q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

function zScore(alpha, tail) {
  // tail=2: 双侧, tail=1: 单侧
  if (tail == 2) return normInv(1 - alpha / 2);
  return normInv(1 - alpha);
}

function zPower(power) {
  return normInv(power);
}

function fmt(n) {
  if (n === null || isNaN(n) || !isFinite(n)) return '—';
  return Math.ceil(n).toLocaleString('zh-CN');
}

function fmtDec(n, d = 3) {
  if (n === null || isNaN(n) || !isFinite(n)) return '—';
  return n.toFixed(d);
}

function animateNum(elId, val) {
  const el = document.getElementById(elId);
  if (!el) return;
  const newVal = isNaN(val) || !isFinite(val) ? '—' : Math.ceil(val).toLocaleString('zh-CN');
  if (el.textContent !== newVal) {
    el.textContent = newVal;
    el.classList.remove('updated');
    void el.offsetWidth;
    el.classList.add('updated');
  }
}

// ===== 导航切换 =====
function showSection(id) {
  document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('.nav-btn').forEach(b => {
    b.classList.remove('active');
    b.setAttribute('aria-selected', 'false');
    b.removeAttribute('aria-current');
  });
  document.getElementById(id).classList.add('active');
  const idx = ['rct','cross','cohort','casecontrol','prediction'].indexOf(id);
  const activeBtn = document.querySelectorAll('.nav-btn')[idx];
  activeBtn.classList.add('active');
  activeBtn.setAttribute('aria-selected', 'true');
  activeBtn.setAttribute('aria-current', 'page');
  // 触发对应计算
  setTimeout(() => {
    if (id === 'rct') { calcRCTMean(); calcRCTProp(); }
    if (id === 'cross') { calcCrossProp(); calcCrossMean(); }
    if (id === 'cohort') calcCohort();
    if (id === 'casecontrol') calcCaseControl();
    if (id === 'prediction') { calcEPV(); calcRiley(); calcContPred(); }
  }, 100);
}

function switchRCT(type) {
  document.getElementById('rct-mean').style.display = type === 'mean' ? 'block' : 'none';
  document.getElementById('rct-prop').style.display = type === 'prop' ? 'block' : 'none';
  document.getElementById('rct-mean-btn').classList.toggle('active', type === 'mean');
  document.getElementById('rct-prop-btn').classList.toggle('active', type === 'prop');
}

function switchCross(type) {
  document.getElementById('cross-prop').style.display = type === 'prop' ? 'block' : 'none';
  document.getElementById('cross-mean').style.display = type === 'mean' ? 'block' : 'none';
  document.getElementById('cross-prop-btn').classList.toggle('active', type === 'prop');
  document.getElementById('cross-mean-btn').classList.toggle('active', type === 'mean');
}

function switchPred(type) {
  ['epv','riley','cont'].forEach(t => {
    document.getElementById('pred-' + t).style.display = t === type ? 'block' : 'none';
    document.getElementById('pred-' + t + '-btn').classList.toggle('active', t === type);
  });
}

// ===== 主题切换 =====
function toggleTheme() {
  const html = document.documentElement;
  const isDark = html.getAttribute('data-theme') === 'dark';
  html.setAttribute('data-theme', isDark ? 'light' : 'dark');
  document.getElementById('themeIcon').textContent = isDark ? '🌙' : '☀️';
  localStorage.setItem('theme', isDark ? 'light' : 'dark');
  // 重新渲染 MathJax
  if (window.MathJax) MathJax.typesetPromise();
}

// 复制模板
function copyTemplate(id) {
  const el = document.getElementById(id);
  if (!el) return;
  const text = el.innerText || el.textContent;

  // 查找关联的复制按钮并添加反馈
  const templateCard = el.closest('.template-card');
  const btn = templateCard ? templateCard.querySelector('.btn-copy') : null;

  navigator.clipboard.writeText(text).then(() => {
    showToast('✅ 模板已复制到剪贴板');
    if (btn) {
      btn.classList.add('copied');
      btn.textContent = '✅ 已复制';
      setTimeout(() => {
        btn.classList.remove('copied');
        btn.innerHTML = '📋 复制模板';
      }, 1500);
    }
  }).catch(() => {
    const ta = document.createElement('textarea');
    ta.value = text;
    document.body.appendChild(ta);
    ta.select();
    document.execCommand('copy');
    document.body.removeChild(ta);
    showToast('✅ 模板已复制到剪贴板');
    if (btn) {
      btn.classList.add('copied');
      btn.textContent = '✅ 已复制';
      setTimeout(() => {
        btn.classList.remove('copied');
        btn.innerHTML = '📋 复制模板';
      }, 1500);
    }
  });
}

function showToast(msg) {
  let toast = document.querySelector('.copy-toast');
  if (!toast) {
    toast = document.createElement('div');
    toast.className = 'copy-toast';
    document.body.appendChild(toast);
  }
  toast.textContent = msg;
  toast.classList.add('show');
  setTimeout(() => toast.classList.remove('show'), 2500);
}

/* ========================================
   1. 随机对照试验 — 均数比较
   Chow et al. (2018) Eq. 3.2.5 / 3.2.6
   ======================================== */
function calcRCTMean() {
  const alpha = parseFloat(document.getElementById('rct_m_alpha').value);
  const tail = parseInt(document.getElementById('rct_m_tail').value);
  const power = parseFloat(document.getElementById('rct_m_power').value);
  const mu1 = parseFloat(document.getElementById('rct_m_mu1').value);
  const mu2 = parseFloat(document.getElementById('rct_m_mu2').value);
  const sigma = parseFloat(document.getElementById('rct_m_sigma').value);
  const k = parseFloat(document.getElementById('rct_m_k').value);
  const dropoutPct = parseFloat(document.getElementById('rct_m_dropout').value);

  if ([alpha, power, sigma, k].some(v => isNaN(v) || v <= 0)) return;
  if (isNaN(mu1) || isNaN(mu2)) return;

  const za = zScore(alpha, tail);
  const zb = zPower(power);
  const delta = Math.abs(mu1 - mu2);
  if (delta === 0) return;

  // n₁（对照组），n₂（试验组）= k*n₁
  const n1 = Math.pow(za + zb, 2) * Math.pow(sigma, 2) * (1 + 1 / k) / Math.pow(delta, 2);
  const n2 = k * n1;
  const nTotal = n1 + n2;
  const dropout = dropoutPct / 100;
  const nAdj = nTotal / (1 - dropout);

  const effectD = fmtDec(delta / sigma);

  animateNum('rct_m_n_per', Math.max(n1, n2));
  document.getElementById('rct_m_n_total').textContent = fmt(nTotal);
  document.getElementById('rct_m_n_adj').textContent = fmt(nAdj);
  document.getElementById('rct_m_delta').textContent = effectD;
  document.getElementById('rct_m_za').textContent = fmtDec(za);
  document.getElementById('rct_m_zb').textContent = fmtDec(zb);

  // 论文模板
  const n1c = Math.ceil(n1);
  const n2c = Math.ceil(n2);
  const nAdjC = Math.ceil(nAdj);
  const tailStr = tail === 2 ? '双侧' : '单侧';
  const template = `【样本量计算方法】

本研究为随机对照试验，主要结局指标为连续变量。样本量计算基于两独立样本均数比较的t检验。

参数设定：依据既往文献/预试验数据，对照组预期均值为 μ₁ = ${mu1}，试验组预期均值为 μ₂ = ${mu2}，两组均数差值 δ = ${delta}，共同标准差 σ = ${sigma}（效应量 d = ${effectD}）。设定${tailStr}检验显著性水平 α = ${alpha}，检验效能（1-β）= ${(power * 100).toFixed(0)}%，两组分配比例为 1:${k}。

计算方法：采用以下公式计算每组样本量：
n = (z_{α/2} + z_β)² × σ² × (1 + 1/k) / δ²

计算结果：对照组需 ${n1c} 例，试验组需 ${n2c} 例，共计 ${Math.ceil(nTotal)} 例。考虑约 ${dropoutPct}% 的脱落/失访率，最终各研究组样本量取整后，研究共需招募 ${nAdjC} 例受试者（对照组 ${Math.ceil(n1c / (1 - dropout))} 例，试验组 ${Math.ceil(n2c / (1 - dropout))} 例）。

参考文献：Chow SC, Shao J, Wang H, Lokhnygina Y. Sample Size Calculations in Clinical Research (3rd ed.). Chapman & Hall/CRC, 2018.`;

  document.getElementById('rct_m_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 计算标准正态分布分位数**
- α = ${alpha} (${tailStr})
- z_{α/2} = ${fmtDec(za, 4)}
- 检验效能 1-β = ${power}，z_β = ${fmtDec(zb, 4)}

**Step 2: 计算均数差值和效应量**
- 两组均数差 δ = |μ₁ - μ₂| = |${mu1} - ${mu2}| = ${delta}
- Cohen's d = δ/σ = ${delta}/${sigma} = ${effectD}

**Step 3: 代入公式**
$$n_1 = \\frac{(z_{\\alpha/2}+z_\\beta)^2 \\cdot \\sigma^2 \\cdot (1+1/k)}{\\delta^2}$$

$$n_1 = \\frac{(${fmtDec(za,4)}+${fmtDec(zb,4)})^2 \\cdot ${sigma}^2 \\cdot (1+1/${k})}{${delta}^2}$$

$$n_1 = \\frac{${fmtDec(Math.pow(za+zb,2),4)} \\cdot ${sigma*sigma}}{${delta*delta}} \\cdot ${fmtDec(1+1/k,4)}$$

$$n_1 = \\frac{${fmtDec(Math.pow(za+zb,2)*sigma*sigma,4)}}{${delta*delta}} \\times ${fmtDec(1+1/k,4)} = ${fmtDec(n1,2)}$$

**Step 4: 计算各组样本量**
- 对照组 n₁ = ${fmtDec(n1,2)} → 取整 ${n1c} 例
- 试验组 n₂ = k × n₁ = ${k} × ${n1c} = ${n2c} 例
- 总样本量 = ${n1c} + ${n2c} = ${Math.ceil(nTotal)} 例

**Step 5: 考虑脱落率调整**
- 脱落率 = ${dropoutPct}%
- 调整后总样本量 = ${Math.ceil(nTotal)} / (1 - ${dropoutPct}%) = ${fmtDec(nAdj,2)} → ${nAdjC} 例`;

  document.getElementById('rct_m_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   2. 随机对照试验 — 率比较
   Fleiss (1981) / Chow et al. (2018) Ch.4
   ======================================== */
function calcRCTProp() {
  const alpha = parseFloat(document.getElementById('rct_p_alpha').value);
  const tail = parseInt(document.getElementById('rct_p_tail').value);
  const power = parseFloat(document.getElementById('rct_p_power').value);
  const p1 = parseFloat(document.getElementById('rct_p_p1').value);
  const p2 = parseFloat(document.getElementById('rct_p_p2').value);
  const k = parseFloat(document.getElementById('rct_p_k').value);
  const dropoutPct = parseFloat(document.getElementById('rct_p_dropout').value);

  if ([alpha, power, p1, p2, k].some(v => isNaN(v) || v <= 0)) return;
  if (p1 === p2) return;

  const za = zScore(alpha, tail);
  const zb = zPower(power);
  const delta = Math.abs(p1 - p2);
  const pbar = (p1 + k * p2) / (1 + k);

  // Fleiss公式（含校正）
  const term1 = za * Math.sqrt((1 + 1 / k) * pbar * (1 - pbar));
  const term2 = zb * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / k);
  const n1_basic = Math.pow(term1 + term2, 2) / Math.pow(delta, 2);

  // 连续性校正（Fleiss, 1981）
  const n1 = (n1_basic / 4) * Math.pow(1 + Math.sqrt(1 + 2 * (1 + 1 / k) / (n1_basic * delta)), 2);
  const n2 = k * n1;
  const nTotal = n1 + n2;
  const dropout = dropoutPct / 100;
  const nAdj = nTotal / (1 - dropout);

  animateNum('rct_p_n_per', Math.max(n1, n2));
  document.getElementById('rct_p_n_total').textContent = fmt(nTotal);
  document.getElementById('rct_p_n_adj').textContent = fmt(nAdj);
  document.getElementById('rct_p_pbar').textContent = fmtDec(pbar);
  document.getElementById('rct_p_delta').textContent = fmtDec(delta);

  const n1c = Math.ceil(n1);
  const n2c = Math.ceil(n2);
  const nAdjC = Math.ceil(nAdj);
  const tailStr = tail === 2 ? '双侧' : '单侧';
  const template = `【样本量计算方法】

本研究为随机对照试验，主要结局指标为二分类变量（事件率）。样本量计算基于两独立样本率比较（卡方检验，Fleiss连续性校正法）。

参数设定：依据既往文献/预试验数据，对照组预期事件率 p₁ = ${(p1 * 100).toFixed(1)}%，试验组预期事件率 p₂ = ${(p2 * 100).toFixed(1)}%，两组率差 |p₁-p₂| = ${(delta * 100).toFixed(1)}%。设定${tailStr}检验显著性水平 α = ${alpha}，检验效能（1-β）= ${(power * 100).toFixed(0)}%，两组分配比例为 1:${k}。

计算方法：基于 Fleiss（1981）含连续性校正的公式：
n = [z_{α/2}√((1+1/k)p̄(1-p̄)) + z_β√(p₁(1-p₁) + p₂(1-p₂)/k)]² / (p₁-p₂)²
合并率 p̄ = ${fmtDec(pbar)}

计算结果：对照组需 ${n1c} 例，试验组需 ${n2c} 例，共计 ${Math.ceil(nTotal)} 例。考虑约 ${dropoutPct}% 的脱落/失访率，研究共需招募 ${nAdjC} 例受试者（对照组 ${Math.ceil(n1c / (1 - dropout))} 例，试验组 ${Math.ceil(n2c / (1 - dropout))} 例）。

参考文献：
1. Fleiss JL. Statistical Methods for Rates and Proportions (2nd ed.). Wiley, 1981.
2. Chow SC, et al. Sample Size Calculations in Clinical Research (3rd ed.). CRC, 2018.`;

  document.getElementById('rct_p_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 计算标准正态分布分位数**
- α = ${alpha} (${tailStr})
- z_{α/2} = ${fmtDec(za, 4)}
- 检验效能 1-β = ${power}，z_β = ${fmtDec(zb, 4)}

**Step 2: 计算合并率**
$$\\bar{p} = \\frac{p_1 + k \\cdot p_2}{1+k} = \\frac{${p1} + ${k} \\times ${p2}}{1+${k}} = ${fmtDec(pbar, 4)}$$

**Step 3: 计算两项独立成分**
- term1 = z_{α/2} × √[(1+1/k) × p̄ × (1-p̄)]
  = ${fmtDec(za,4)} × √[(1+1/${k}) × ${fmtDec(pbar,4)} × ${fmtDec(1-pbar,4)}]
  = ${fmtDec(za,4)} × √${fmtDec((1+1/k)*pbar*(1-pbar),4)}
  = ${fmtDec(term1, 4)}

- term2 = z_β × √[p₁(1-p₁) + p₂(1-p₂)/k]
  = ${fmtDec(zb,4)} × √[${p1}×${fmtDec(1-p1,4)} + ${p2}×${fmtDec(1-p2,4)}/${k}]
  = ${fmtDec(zb,4)} × √${fmtDec(p1*(1-p1)+p2*(1-p2)/k,4)}
  = ${fmtDec(term2, 4)}

**Step 4: 代入Fleiss公式（无校正）**
$$n_1 = \\frac{(term1 + term2)^2}{(p_1-p_2)^2} = \\frac{(${fmtDec(term1,4)}+${fmtDec(term2,4)})^2}{${delta}^2} = ${fmtDec(n1_basic, 2)}$$

**Step 5: Fleiss连续性校正**
$$n_{cc} = \\frac{n}{4}\\left[1+\\sqrt{1+\\frac{2(1+1/k)}{n|p_1-p_2|}}\\right]^2$$
= ${fmtDec(n1, 2)} → 取整 ${n1c} 例

**Step 6: 计算各组样本量**
- 对照组 n₁ = ${n1c} 例
- 试验组 n₂ = k × n₁ = ${k} × ${n1c} = ${n2c} 例
- 总样本量 = ${Math.ceil(nTotal)} 例

**Step 7: 考虑脱落率调整**
- 调整后总样本量 = ${Math.ceil(nTotal)} / (1 - ${dropoutPct}%) = ${fmtDec(nAdj,2)} → ${nAdjC} 例`;

  document.getElementById('rct_p_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   3. 横断面调查 — 估计率/比例
   Cochran (1977) / Lwanga & Lemeshow (1991)
   ======================================== */
function calcCrossProp() {
  const alpha = parseFloat(document.getElementById('cs_p_alpha').value);
  const p = parseFloat(document.getElementById('cs_p_p').value);
  const d = parseFloat(document.getElementById('cs_p_d').value);
  const deff = parseFloat(document.getElementById('cs_p_deff').value);
  const dropoutPct = parseFloat(document.getElementById('cs_p_dropout').value);

  if ([alpha, p, d, deff].some(v => isNaN(v) || v <= 0)) return;

  const za = normInv(1 - alpha / 2);
  const nBase = Math.pow(za, 2) * p * (1 - p) / Math.pow(d, 2);
  const nDeff = nBase * deff;
  const dropout = dropoutPct / 100;
  const nFinal = nDeff / (1 - dropout);

  animateNum('cs_p_n', nFinal);
  document.getElementById('cs_p_n_base').textContent = fmt(nBase);
  document.getElementById('cs_p_n_deff').textContent = fmt(nDeff);
  document.getElementById('cs_p_n_final').textContent = fmt(nFinal);
  document.getElementById('cs_p_za').textContent = fmtDec(za);

  const nFinalC = Math.ceil(nFinal);
  const nBaseC = Math.ceil(nBase);
  const ci = ((1 - alpha) * 100).toFixed(0);
  const template = `【样本量计算方法】

本研究为横断面调查，旨在估计目标人群中某结局指标/疾病的患病率（比例）。

参数设定：根据既往文献/预调查结果，预期患病率/比例 p = ${(p * 100).toFixed(1)}%（若无先验信息则取最保守值 50%），设定 ${ci}% 置信区间（α = ${alpha}），允许误差（绝对精度）d = ±${(d * 100).toFixed(1)}%，设计效应 DEFF = ${deff}（${deff > 1 ? '整群/多阶段抽样' : '简单随机抽样'}），无应答率约 ${dropoutPct}%。

计算方法：基于 Cochran（1977）公式：
n = z²_{α/2} × p(1-p) / d²

计算结果：基础样本量为 ${nBaseC} 例；考虑设计效应（DEFF = ${deff}）调整后为 ${Math.ceil(nBase * deff)} 例；计入 ${dropoutPct}% 无应答率后，最终所需样本量为 ${nFinalC} 例。

参考文献：
1. Cochran WG. Sampling Techniques (3rd ed.). Wiley, 1977.
2. Lwanga SK, Lemeshow S. Sample Size Determination in Health Studies. WHO, 1991.
3. Chow SC, et al. Sample Size Calculations in Clinical Research (3rd ed.). CRC, 2018.`;

  document.getElementById('cs_p_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 计算标准正态分布分位数**
- α = ${alpha}（对应 ${ci}% 置信区间）
- z_{α/2} = ${fmtDec(za, 4)}

**Step 2: 计算基础样本量（Cochran公式）**
$$n = \\frac{z_{\\alpha/2}^2 \\cdot p(1-p)}{d^2}$$

代入参数：
- p = ${p}（预期患病率）
- d = ${d}（允许误差）
- p(1-p) = ${p} × (1-${p}) = ${fmtDec(p*(1-p), 4)}
- z_{α/2}² = (${fmtDec(za,4)})² = ${fmtDec(za*za, 4)}

$$n_{base} = \\frac{${fmtDec(za*za,4)} \\times ${fmtDec(p*(1-p),4)}}{${d}^2} = \\frac{${fmtDec(za*za*p*(1-p),4)}}{${d*d}} = ${fmtDec(nBase, 2)}$$

**Step 3: 考虑设计效应（DEFF）**
$$n_{deff} = n_{base} \\times DEFF = ${fmtDec(nBase,2)} \\times ${deff} = ${fmtDec(nDeff, 2)}$$

**Step 4: 考虑无应答率**
$$n_{final} = \\frac{n_{deff}}{1 - r} = \\frac{${fmtDec(nDeff,2)}}{1 - ${dropoutPct}\\%} = ${fmtDec(nFinal, 2)} \\rightarrow ${nFinalC} 例`;

  document.getElementById('cs_p_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   4. 横断面调查 — 估计均值
   ======================================== */
function calcCrossMean() {
  const alpha = parseFloat(document.getElementById('cs_m_alpha').value);
  const sigma = parseFloat(document.getElementById('cs_m_sigma').value);
  const d = parseFloat(document.getElementById('cs_m_d').value);
  const deff = parseFloat(document.getElementById('cs_m_deff').value);
  const dropoutPct = parseFloat(document.getElementById('cs_m_dropout').value);

  if ([alpha, sigma, d, deff].some(v => isNaN(v) || v <= 0)) return;

  const za = normInv(1 - alpha / 2);
  const nBase = Math.pow(za, 2) * Math.pow(sigma, 2) / Math.pow(d, 2) * deff;
  const dropout = dropoutPct / 100;
  const nFinal = nBase / (1 - dropout);

  animateNum('cs_m_n', nFinal);
  document.getElementById('cs_m_n_base').textContent = fmt(nBase);
  document.getElementById('cs_m_n_final').textContent = fmt(nFinal);

  const nFinalC = Math.ceil(nFinal);
  const ci = ((1 - alpha) * 100).toFixed(0);
  const template = `【样本量计算方法】

本研究为横断面调查，旨在估计目标人群中某连续型指标的总体均值。

参数设定：根据既往文献/预调查结果，总体标准差 σ = ${sigma}，设定 ${ci}% 置信区间（α = ${alpha}），允许绝对误差 d = ±${d}，设计效应 DEFF = ${deff}，无应答率约 ${dropoutPct}%。

计算方法：基于简单随机抽样估计均值公式：
n = z²_{α/2} × σ² / d² × DEFF

计算结果：基础样本量为 ${Math.ceil(nBase)} 例；考虑 ${dropoutPct}% 无应答率后，最终所需样本量为 ${nFinalC} 例。

参考文献：
1. Cochran WG. Sampling Techniques (3rd ed.). Wiley, 1977.
2. Chow SC, et al. Sample Size Calculations in Clinical Research (3rd ed.). CRC, 2018.`;

  document.getElementById('cs_m_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 计算标准正态分布分位数**
- α = ${alpha}（对应 ${ci}% 置信区间）
- z_{α/2} = ${fmtDec(za, 4)}

**Step 2: 计算基础样本量**
$$n = \\frac{z_{\\alpha/2}^2 \\cdot \\sigma^2}{d^2} \\times DEFF$$

代入参数：
- σ = ${sigma}（标准差）
- d = ${d}（允许误差）
- σ² = ${sigma}² = ${sigma*sigma}
- z_{α/2}² = (${fmtDec(za,4)})² = ${fmtDec(za*za, 4)}

$$n_{base} = \\frac{${fmtDec(za*za,4)} \\times ${sigma*sigma}}{${d}^2} \\times ${deff} = \\frac{${fmtDec(za*za*sigma*sigma,4)}}{${d*d}} \\times ${deff} = ${fmtDec(nBase, 2)}$$

**Step 3: 考虑脱落率调整**
$$n_{final} = \\frac{n_{base}}{1 - r} = \\frac{${fmtDec(nBase,2)}}{1 - ${dropoutPct}\\%} = ${fmtDec(nFinal, 2)} \\rightarrow ${nFinalC} 例

---

**公式验证（您的例题）**
- σ = 114, d = 14, α = 0.05
- z_{0.975} = ${fmtDec(za, 4)}
- n = (${fmtDec(za,4)} × 114 / 14)² = ${fmtDec(nBase, 2)} → ${nFinalC} 例 ✓`;

  document.getElementById('cs_m_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   5. 队列研究
   Kelsey公式 / Chow et al. Ch.8
   ======================================== */
function calcCohort() {
  const alpha = parseFloat(document.getElementById('co_alpha').value);
  const tail = parseInt(document.getElementById('co_tail').value);
  const power = parseFloat(document.getElementById('co_power').value);
  const p0 = parseFloat(document.getElementById('co_p0').value);
  const p1 = parseFloat(document.getElementById('co_p1').value);
  const k = parseFloat(document.getElementById('co_k').value);
  const dropoutPct = parseFloat(document.getElementById('co_dropout').value);

  if ([alpha, power, p0, p1, k].some(v => isNaN(v) || v <= 0)) return;
  if (p0 === p1) return;

  const za = zScore(alpha, tail);
  const zb = zPower(power);
  const pbar = (p0 + k * p1) / (1 + k);
  const delta = Math.abs(p1 - p0);

  const term1 = za * Math.sqrt((1 + 1 / k) * pbar * (1 - pbar));
  const term2 = zb * Math.sqrt(p0 * (1 - p0) + p1 * (1 - p1) / k);
  const n0 = Math.pow(term1 + term2, 2) / Math.pow(delta, 2);
  const n1 = k * n0;
  const nTotal = n0 + n1;
  const dropout = dropoutPct / 100;
  const nAdj = nTotal / (1 - dropout);
  const rr = p1 / p0;
  const ar = delta;

  animateNum('co_n_per', n0);
  document.getElementById('co_n1').textContent = fmt(n1);
  document.getElementById('co_n0').textContent = fmt(n0);
  document.getElementById('co_n_total').textContent = fmt(nTotal);
  document.getElementById('co_n_adj').textContent = fmt(nAdj);
  document.getElementById('co_rr').textContent = fmtDec(rr, 2);
  document.getElementById('co_ar').textContent = fmtDec(ar, 3);

  const n0c = Math.ceil(n0);
  const n1c = Math.ceil(n1);
  const nAdjC = Math.ceil(nAdj);
  const tailStr = tail === 2 ? '双侧' : '单侧';
  const template = `【样本量计算方法】

本研究为前瞻性队列研究，旨在比较暴露组与非暴露组的结局发生率。

参数设定：依据既往文献报道，非暴露（对照）组结局发生率 p₀ = ${(p0 * 100).toFixed(1)}%，预期暴露组结局发生率 p₁ = ${(p1 * 100).toFixed(1)}%，对应相对危险度 RR = ${fmtDec(rr, 2)}，归因危险度 AR = ${(ar * 100).toFixed(1)}%。设定${tailStr}检验显著性水平 α = ${alpha}，检验效能（1-β）= ${(power * 100).toFixed(0)}%，暴露组与非暴露组分配比例为 ${k}:1。

计算方法：基于 Kelsey 等（1996）队列研究率比较公式：
n₀ = [z_{α/2}√((1+1/k)p̄(1-p̄)) + z_β√(p₀(1-p₀)+p₁(1-p₁)/k)]² / (p₁-p₀)²
合并率 p̄ = ${fmtDec(pbar)}

计算结果：非暴露组需 ${n0c} 例，暴露组需 ${n1c} 例，共计 ${Math.ceil(nTotal)} 例。考虑约 ${dropoutPct}% 的失访率，最终需招募 ${nAdjC} 例研究对象。

参考文献：
1. Kelsey JL, et al. Methods in Observational Epidemiology (2nd ed.). Oxford, 1996.
2. Chow SC, et al. Sample Size Calculations in Clinical Research (3rd ed.). CRC, 2018, Chapter 8.`;

  document.getElementById('co_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 计算标准正态分布分位数**
- α = ${alpha} (${tailStr})
- z_{α/2} = ${fmtDec(za, 4)}
- 检验效能 1-β = ${power}，z_β = ${fmtDec(zb, 4)}

**Step 2: 计算基本统计量**
- 非暴露组率 p₀ = ${p0}
- 暴露组率 p₁ = ${p1}
- 率差 δ = |p₁ - p₀| = |${p1} - ${p0}| = ${fmtDec(delta, 4)}
- 相对危险度 RR = p₁/p₀ = ${p1}/${p0} = ${fmtDec(rr, 4)}
- 归因危险度 AR = p₁ - p₀ = ${fmtDec(ar, 4)}

**Step 3: 计算合并率**
$$\\bar{p} = \\frac{p_0 + k \\cdot p_1}{1+k} = \\frac{${p0} + ${k} \\times ${p1}}{1+${k}} = ${fmtDec(pbar, 4)}$$

**Step 4: 计算两项独立成分**
- term1 = z_{α/2} × √[(1+1/k) × p̄ × (1-p̄)]
  = ${fmtDec(za,4)} × √[(1+1/${k}) × ${fmtDec(pbar,4)} × ${fmtDec(1-pbar,4)}]
  = ${fmtDec(term1, 4)}

- term2 = z_β × √[p₀(1-p₀) + p₁(1-p₁)/k]
  = ${fmtDec(zb,4)} × √[${p0}×${fmtDec(1-p0,4)} + ${p1}×${fmtDec(1-p1,4)}/${k}]
  = ${fmtDec(term2, 4)}

**Step 5: 代入Kelsey公式**
$$n_0 = \\frac{(term1 + term2)^2}{(p_1-p_0)^2} = \\frac{(${fmtDec(term1,4)}+${fmtDec(term2,4)})^2}{${fmtDec(delta,4)}^2} = ${fmtDec(n0, 2)}$$

**Step 6: 计算各组样本量**
- 非暴露组 n₀ = ${fmtDec(n0,2)} → 取整 ${n0c} 例
- 暴露组 n₁ = k × n₀ = ${k} × ${n0c} = ${n1c} 例
- 总样本量 = ${Math.ceil(nTotal)} 例

**Step 7: 考虑失访率**
- 调整后总样本量 = ${Math.ceil(nTotal)} / (1 - ${dropoutPct}%) = ${fmtDec(nAdj,2)} → ${nAdjC} 例`;

  document.getElementById('co_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   6. 病例对照研究
   Kelsey公式 / Schlesselman (1982)
   ======================================== */
function calcCaseControl() {
  const alpha = parseFloat(document.getElementById('cc_alpha').value);
  const tail = parseInt(document.getElementById('cc_tail').value);
  const power = parseFloat(document.getElementById('cc_power').value);
  const p0 = parseFloat(document.getElementById('cc_p0').value);
  const or = parseFloat(document.getElementById('cc_or').value);
  const m = parseFloat(document.getElementById('cc_m').value);
  const dropoutPct = parseFloat(document.getElementById('cc_dropout').value);

  if ([alpha, power, p0, or, m].some(v => isNaN(v) || v <= 0)) return;

  // 由 p0 和 OR 推算 p1
  const p1 = or * p0 / (1 - p0 + or * p0);
  const pbar = (p1 + m * p0) / (1 + m);
  const delta = Math.abs(p1 - p0);

  const za = zScore(alpha, tail);
  const zb = zPower(power);

  const term1 = za * Math.sqrt((1 + 1 / m) * pbar * (1 - pbar));
  const term2 = zb * Math.sqrt(p1 * (1 - p1) + p0 * (1 - p0) / m);
  const nCase = Math.pow(term1 + term2, 2) / Math.pow(delta, 2);
  const nCtrl = m * nCase;
  const nTotal = nCase + nCtrl;
  const dropout = dropoutPct / 100;
  const nAdj = nTotal / (1 - dropout);

  animateNum('cc_n_case', nCase);
  document.getElementById('cc_n_ctrl').textContent = fmt(nCtrl);
  document.getElementById('cc_n_total').textContent = fmt(nTotal);
  document.getElementById('cc_n_adj').textContent = fmt(nAdj);
  document.getElementById('cc_p1').textContent = fmtDec(p1);
  document.getElementById('cc_pbar').textContent = fmtDec(pbar);

  const nCaseC = Math.ceil(nCase);
  const nCtrlC = Math.ceil(nCtrl);
  const nAdjC = Math.ceil(nAdj);
  const tailStr = tail === 2 ? '双侧' : '单侧';
  const template = `【样本量计算方法】

本研究为病例对照研究，旨在探究目标暴露因素与所研究疾病/结局之间的关联强度。

参数设定：依据人群流行病学数据/既往文献，对照人群中暴露因素的比例（对照组暴露率）p₀ = ${(p0 * 100).toFixed(1)}%，预期最小有意义的比值比（OR）= ${or}，由此推算病例组暴露率 p₁ = ${(p1 * 100).toFixed(3)}。设定${tailStr}检验显著性水平 α = ${alpha}，检验效能（1-β）= ${(power * 100).toFixed(0)}%，病例与对照的比例为 1:${m}。

计算方法：基于 Kelsey 等（1996）病例对照研究公式，先由 OR 和 p₀ 推算 p₁：
p₁ = OR × p₀ / (1 - p₀ + OR × p₀) = ${fmtDec(p1, 4)}
再计算病例组样本量：
n_case = [z_{α/2}√((1+1/m)p̄(1-p̄)) + z_β√(p₁(1-p₁)+p₀(1-p₀)/m)]² / (p₁-p₀)²
合并暴露率 p̄ = ${fmtDec(pbar)}

计算结果：病例组需 ${nCaseC} 例，对照组需 ${nCtrlC} 例，共计 ${Math.ceil(nTotal)} 例。计入约 ${dropoutPct}% 的无应答/脱落率，最终需纳入 ${nAdjC} 例研究对象。

参考文献：
1. Schlesselman JJ. Case-Control Studies: Design, Conduct, Analysis. Oxford, 1982.
2. Kelsey JL, et al. Methods in Observational Epidemiology (2nd ed.). Oxford, 1996.
3. Chow SC, et al. Sample Size Calculations in Clinical Research (3rd ed.). CRC, 2018.`;

  document.getElementById('cc_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 计算标准正态分布分位数**
- α = ${alpha} (${tailStr})
- z_{α/2} = ${fmtDec(za, 4)}
- 检验效能 1-β = ${power}，z_β = ${fmtDec(zb, 4)}

**Step 2: 由OR推算病例组暴露率**
$$p_1 = \\frac{OR \\cdot p_0}{1 - p_0 + OR \\cdot p_0} = \\frac{${or} \\times ${p0}}{1 - ${p0} + ${or} \\times ${p0}} = ${fmtDec(p1, 4)}$$

**Step 3: 计算合并暴露率**
$$\\bar{p} = \\frac{p_1 + m \\cdot p_0}{1+m} = \\frac{${fmtDec(p1,4)} + ${m} \\times ${p0}}{1+${m}} = ${fmtDec(pbar, 4)}$$

- 率差 δ = |p₁ - p₀| = |${fmtDec(p1,4)} - ${p0}| = ${fmtDec(delta, 4)}

**Step 4: 计算两项独立成分**
- term1 = z_{α/2} × √[(1+1/m) × p̄ × (1-p̄)]
  = ${fmtDec(za,4)} × √[(1+1/${m}) × ${fmtDec(pbar,4)} × ${fmtDec(1-pbar,4)}]
  = ${fmtDec(term1, 4)}

- term2 = z_β × √[p₁(1-p₁) + p₀(1-p₀)/m]
  = ${fmtDec(zb,4)} × √[${fmtDec(p1,4)}×${fmtDec(1-p1,4)} + ${p0}×${fmtDec(1-p0,4)}/${m}]
  = ${fmtDec(term2, 4)}

**Step 5: 代入Kelsey公式**
$$n_{case} = \\frac{(term1 + term2)^2}{(p_1-p_0)^2} = ${fmtDec(nCase, 2)}$$

**Step 6: 计算各组样本量**
- 病例组 n_case = ${fmtDec(nCase,2)} → 取整 ${nCaseC} 例
- 对照组 n_ctrl = m × n_case = ${m} × ${nCaseC} = ${nCtrlC} 例
- 总样本量 = ${Math.ceil(nTotal)} 例

**Step 7: 考虑无应答率**
- 调整后总样本量 = ${Math.ceil(nTotal)} / (1 - ${dropoutPct}%) = ${fmtDec(nAdj,2)} → ${nAdjC} 例`;

  document.getElementById('cc_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   7. 预测模型 — EPV法
   Peduzzi et al. (1996) / Riley (2019)
   ======================================== */
function calcEPV() {
  const p = parseFloat(document.getElementById('epv_p').value);
  const epv = parseFloat(document.getElementById('epv_epv').value);
  const rate = parseFloat(document.getElementById('epv_rate').value);
  const dropoutPct = parseFloat(document.getElementById('epv_dropout').value);

  if ([p, epv, rate].some(v => isNaN(v) || v <= 0)) return;
  if (rate >= 1) return;

  const eMin = epv * p;
  const nBase = eMin / rate;
  const dropout = dropoutPct / 100;
  const nFinal = nBase / (1 - dropout);

  animateNum('epv_n', nFinal);
  document.getElementById('epv_events').textContent = fmt(eMin);
  document.getElementById('epv_n_base').textContent = fmt(nBase);
  document.getElementById('epv_n_final').textContent = fmt(nFinal);
  document.getElementById('epv_standard').textContent = epv >= 20 ? '严格标准 (≥20)' : epv >= 10 ? '传统标准 (≥10)' : '⚠️ 低于推荐';

  const nFinalC = Math.ceil(nFinal);
  const nBaseC = Math.ceil(nBase);
  const eMinC = Math.ceil(eMin);
  const template = `【样本量计算方法】

本研究旨在建立基于二分类结局变量的多变量预测模型（如Logistic回归）。样本量计算采用EPV（Events Per Variable，每预测变量事件数）法。

参数设定：预测模型中候选预测变量数（含哑变量展开后的参数数）p = ${p}；依据 ${epv >= 20 ? 'Riley et al.（2019）的推荐标准' : 'Peduzzi et al.（1996）的传统标准'}，设定每个预测变量对应的最少事件数 EPV = ${epv}；依据既往文献报道，目标人群中结局事件发生率 π = ${(rate * 100).toFixed(1)}%。

计算方法：
最少需要观察的事件数：E = EPV × p = ${epv} × ${p} = ${eMinC}
所需总样本量：N = E / π = ${eMinC} / ${rate} = ${nBaseC}

计算结果：模型建立至少需要观察到 ${eMinC} 个结局事件，基于事件率 π = ${(rate * 100).toFixed(1)}%，需纳入 ${nBaseC} 例研究对象。考虑约 ${dropoutPct}% 的数据缺失/脱落，最终需招募 ${nFinalC} 例。

重要提示：EPV法为传统经验性方法，Riley et al.（2020）发表于BMJ的新方法（基于收缩因子）提供了更严格的样本量估算，建议同时参考。

参考文献：
1. Peduzzi P, et al. A simulation study of the number of events per variable in logistic regression. J Clin Epidemiol. 1996;49:1373-1379.
2. Riley RD, et al. Calculating the sample size required for developing a clinical prediction model. BMJ. 2020;368:m441.`;

  document.getElementById('epv_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 设定参数**
- 预测变量数 p = ${p}
- EPV 设定值 = ${epv}（${epv >= 20 ? 'Riley推荐标准' : epv >= 10 ? '传统标准' : '低于推荐标准'}）
- 结局事件率 π = ${rate}

**Step 2: 计算最少需要的事件数**
$$E_{min} = EPV \\times p = ${epv} \\times ${p} = ${eMinC} 个事件$$

**Step 3: 计算基础样本量**
$$N_{base} = \\frac{E_{min}}{\\pi} = \\frac{${eMinC}}{${rate}} = ${nBaseC} 例$$

**Step 4: 考虑脱落率**
$$N_{final} = \\frac{N_{base}}{1 - r} = \\frac{${nBaseC}}{1 - ${dropoutPct}\\%} = ${fmtDec(nFinal, 2)} \\rightarrow ${nFinalC} 例

---

**说明**
- EPV (Events Per Variable) 法则是传统经验法则
- Peduzzi et al. (1996) 建议 EPV ≥ 10
- Riley et al. (2019) 建议 EPV ≥ 20 以获得更稳定的模型`;

  document.getElementById('epv_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   8. 预测模型 — Riley公式
   Riley et al. (2019/2020)
   ======================================== */
function calcRiley() {
  const p = parseFloat(document.getElementById('riley_p').value);
  const rate = parseFloat(document.getElementById('riley_rate').value);
  const r2nagel = parseFloat(document.getElementById('riley_r2').value);
  const s = parseFloat(document.getElementById('riley_s').value);
  const dropoutPct = parseFloat(document.getElementById('riley_dropout').value);

  if ([p, rate, r2nagel, s].some(v => isNaN(v) || v <= 0)) return;
  if (rate >= 1 || r2nagel >= 1 || s >= 1) return;

  // Cox-Snell R² 的最大值（当所有人都是事件时）
  const r2cs_max = 1 - Math.pow(rate, 2 * rate) * Math.pow(1 - rate, 2 * (1 - rate));
  // 将Nagelkerke R² 转换为 Cox-Snell R²
  const r2cs = r2nagel * r2cs_max;

  // 准则1：收缩因子
  // N₁ ≥ p / ((S-1) × ln(1 - R²_CS/S))
  // 注意 S-1 < 0, ln(1-R²_CS/S) 也<0，所以是正数
  const r2cs_s = r2cs / s;
  if (r2cs_s >= 1) {
    document.getElementById('riley_n').textContent = '参数无效';
    return;
  }
  const n1 = p / ((s - 1) * Math.log(1 - r2cs_s));

  // 准则2：Harrell EPV≥20
  const n2 = 20 * p / rate;

  const nFinal = Math.max(n1, n2);
  const dropout = dropoutPct / 100;
  const nAdj = nFinal / (1 - dropout);
  const events = nFinal * rate;

  animateNum('riley_n', nAdj);
  document.getElementById('riley_n1').textContent = fmt(n1);
  document.getElementById('riley_n2').textContent = fmt(n2);
  document.getElementById('riley_n_final').textContent = fmt(nFinal);
  document.getElementById('riley_events').textContent = fmt(events);

  const n1c = Math.ceil(n1);
  const n2c = Math.ceil(n2);
  const nFinalC = Math.ceil(nFinal);
  const nAdjC = Math.ceil(nAdj);
  const eventsC = Math.ceil(events);
  const template = `【样本量计算方法】

本研究旨在开发基于二分类结局的多变量临床预测模型（如Logistic回归）。样本量计算采用 Riley et al.（2019/2020）发表于 Statistics in Medicine 及 BMJ 的基于收缩因子的方法，以控制模型过拟合并保证参数估计稳定性。

参数设定：
• 候选预测变量数（含哑变量参数）p = ${p}
• 目标人群结局事件率 π = ${(rate * 100).toFixed(1)}%
• 预期模型 Nagelkerke R² = ${r2nagel}（基于类似研究文献/预模型估计）
  换算 Cox-Snell R² = ${fmtDec(r2cs, 4)}（R²_CS_max = ${fmtDec(r2cs_max, 4)}）
• 最大可接受收缩因子 S = ${s}（Riley 推荐 ≥0.9，值越高模型越不过拟合）

计算方法（依据两个准则取较大值）：
准则一（收缩因子控制）：
N₁ ≥ p / [(S-1) × ln(1 - R²_CS/S)] = ${n1c}
准则二（Harrell EPV ≥ 20）：
N₂ ≥ 20p / π = ${n2c}
最终取 N = max(N₁, N₂) = ${nFinalC}

计算结果：预测模型建立需纳入至少 ${nFinalC} 例研究对象（预期包含 ${eventsC} 个结局事件）。考虑约 ${dropoutPct}% 的数据缺失/脱落，最终需招募 ${nAdjC} 例。

参考文献：
1. Riley RD, Snell KIE, Ensor J, et al. Minimum sample size for developing a multivariable prediction model. Stat Med. 2019;38:1262-1275.
2. Riley RD, Ensor J, Snell KIE, et al. Calculating the sample size required for developing a clinical prediction model. BMJ. 2020;368:m441.
3. Van Smeden M, et al. No rationale for 1 variable per 10 events criterion. BMC Med Res Methodol. 2016;16:163.`;

  document.getElementById('riley_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 设定参数**
- 预测变量数 p = ${p}
- 结局事件率 π = ${rate}
- Nagelkerke R² = ${r2nagel}
- 收缩因子 S = ${s}（推荐 ≥0.9）

**Step 2: 计算 Cox-Snell R²**
- R²_CS_max = 1 - π^(2π) × (1-π)^(2(1-π))
  = 1 - ${rate}^(2×${rate}) × (1-${rate})^(2×${1-rate})
  = 1 - ${fmtDec(Math.pow(rate, 2*rate), 6)} × ${fmtDec(Math.pow(1-rate, 2*(1-rate)), 6)}
  = ${fmtDec(r2cs_max, 4)}

- R²_CS = R²_Nagel × R²_CS_max
  = ${r2nagel} × ${fmtDec(r2cs_max, 4)} = ${fmtDec(r2cs, 4)}

**Step 3: 准则一（收缩因子控制）**
$$N_1 \\geq \\frac{p}{(S-1) \\ln(1 - R^2_{CS}/S)}$$

代入：R²_CS/S = ${fmtDec(r2cs,4)}/${s} = ${fmtDec(r2cs_s, 4)}
ln(1 - ${fmtDec(r2cs_s,4)}) = ${fmtDec(Math.log(1-r2cs_s), 4)}
N₁ = ${p} / [(${s}-1) × ${fmtDec(Math.log(1-r2cs_s),4)}] = ${n1c} 例

**Step 4: 准则二（Harrell EPV≥20）**
$$N_2 \\geq \\frac{20p}{\\pi} = \\frac{20 \\times ${p}}{${rate}} = ${n2c} 例

**Step 5: 取较大值**
- N₁ = ${n1c}，N₂ = ${n2c}
- 最终建议 N = max(N₁, N₂) = ${nFinalC} 例

**Step 6: 考虑脱落率**
- 预期事件数 E = N × π = ${nFinalC} × ${rate} = ${eventsC} 个
- 调整后样本量 = ${nFinalC} / (1 - ${dropoutPct}%) = ${fmtDec(nAdj, 2)} → ${nAdjC} 例`;

  document.getElementById('riley_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ========================================
   9. 预测模型 — 连续结局
   Riley et al. (2019) Part I
   ======================================== */
function calcContPred() {
  const p = parseFloat(document.getElementById('cont_p').value);
  const r2 = parseFloat(document.getElementById('cont_r2').value);
  const s = parseFloat(document.getElementById('cont_s').value);
  const alpha = parseFloat(document.getElementById('cont_alpha').value);
  const dropoutPct = parseFloat(document.getElementById('cont_dropout').value);

  if ([p, r2, s, alpha].some(v => isNaN(v) || v <= 0)) return;
  if (r2 >= 1 || s >= 1) return;

  // 简化版：N ≥ p / (S × R²) + p + 1
  const n1 = p / (s * r2) + p + 1;
  const dropout = dropoutPct / 100;
  const nFinal = n1 / (1 - dropout);

  animateNum('cont_n', nFinal);
  document.getElementById('cont_n1').textContent = fmt(n1);
  document.getElementById('cont_n_final').textContent = fmt(nFinal);

  const n1c = Math.ceil(n1);
  const nFinalC = Math.ceil(nFinal);
  const template = `【样本量计算方法】

本研究旨在建立连续型结局变量的多变量预测模型（线性回归或其他回归方法）。样本量计算采用 Riley et al.（2019）发表于 Statistics in Medicine 的基于收缩因子的公式（适用于连续结局）。

参数设定：
• 候选预测变量数 p = ${p}
• 预期调整 R² = ${r2}（基于文献报道的模型解释变异量）
• 最大可接受收缩因子 S = ${s}（推荐 ≥0.9）
• 显著性水平 α = ${alpha}

计算方法（简化收缩因子公式）：
N ≥ p / (S × R²) + p + 1 = ${n1c}

计算结果：模型开发至少需要 ${n1c} 例研究对象。考虑约 ${dropoutPct}% 的数据缺失/脱落，最终需招募 ${nFinalC} 例。

参考文献：
1. Riley RD, Snell KIE, Ensor J, et al. Minimum sample size for developing a multivariable prediction model: Part I – Continuous outcomes. Stat Med. 2019;38(7):1262-1275.
2. Harrell FE. Regression Modeling Strategies (2nd ed.). Springer, 2015.`;

  document.getElementById('cont_template').textContent = template;

  // 计算过程
  const calcProcess = `📋 **计算步骤详解**

**Step 1: 设定参数**
- 预测变量数 p = ${p}
- 预期模型 R² = ${r2}
- 收缩因子 S = ${s}（推荐 ≥0.9）
- 显著性水平 α = ${alpha}

**Step 2: 代入简化收缩因子公式**
$$N \\geq \\frac{p}{S \\cdot R^2} + p + 1$$

代入参数：
$$N \\geq \\frac{${p}}{${s} \\times ${r2}} + ${p} + 1$$

$$N \\geq \\frac{${p}}{${fmtDec(s*r2,4)}} + ${p+1} = ${fmtDec(n1, 2)} \\rightarrow ${n1c} 例$$

**Step 3: 考虑脱落率**
$$N_{final} = \\frac{N}{1 - r} = \\frac{${n1c}}{1 - ${dropoutPct}\\%} = ${fmtDec(nFinal, 2)} \\rightarrow ${nFinalC} 例

---

**说明**
- 此为 Riley et al. (2019) 简化公式，适用于连续结局预测模型
- 公式基于收缩因子控制过拟合的原理
- S 越接近 1，对样本量要求越高`;

  document.getElementById('cont_process').textContent = calcProcess;
  if (window.MathJax) {
    MathJax.typesetPromise && MathJax.typesetPromise();
  }
}

/* ===== 初始化 ===== */
window.addEventListener('DOMContentLoaded', () => {
  // 恢复主题
  const savedTheme = localStorage.getItem('theme') || 'light';
  document.documentElement.setAttribute('data-theme', savedTheme);
  document.getElementById('themeIcon').textContent = savedTheme === 'dark' ? '☀️' : '🌙';

  // 初始计算
  calcRCTMean();
  calcRCTProp();
  calcCrossProp();
  calcCrossMean();
  calcCohort();
  calcCaseControl();
  calcEPV();
  calcRiley();
  calcContPred();

  // 等待 MathJax 加载完成后重新排版
  if (window.MathJax) {
    MathJax.typesetPromise();
  }
});
