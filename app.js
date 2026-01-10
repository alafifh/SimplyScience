// Button ripple
document.addEventListener("click", (e) => {
  const btn = e.target.closest(".btn");
  if (!btn) return;

  const r = document.createElement("span");
  r.className = "ripple";

  const rect = btn.getBoundingClientRect();
  const size = Math.max(rect.width, rect.height);
  r.style.width = r.style.height = size + "px";
  r.style.left = (e.clientX - rect.left - size / 2) + "px";
  r.style.top = (e.clientY - rect.top - size / 2) + "px";

  btn.appendChild(r);
  r.addEventListener("animationend", () => r.remove());
});

// Smooth page-out animation on internal nav (links/buttons that have data-nav)
document.addEventListener("click", (e) => {
  const nav = e.target.closest("[data-nav]");
  if (!nav) return;

  e.preventDefault();
  const href = nav.getAttribute("href") || nav.getAttribute("data-href");
  if (!href) return;

  const page = document.querySelector(".page");
  if (page) page.classList.add("pageOut");

  setTimeout(() => {
    window.location.href = href;
  }, 220);
});

// Simple form “demo” submit (no backend yet)
document.addEventListener("submit", (e) => {
  const form = e.target.closest("form[data-demo]");
  if (!form) return;

  e.preventDefault();
  const msg = form.querySelector("[data-msg]");
  if (msg) {
    msg.textContent = "✅ Submitted (demo). Hook this to your backend later!";
  }
});

function ssGetProfile(){
  try { return JSON.parse(localStorage.getItem("ss_profile") || "null"); }
  catch { return null; }
}

function ssSetProfile(profile){
  localStorage.setItem("ss_profile", JSON.stringify(profile));
}

function ssGetSelectedInterests(){
  const profile = ssGetProfile();
  if (!profile || !Array.isArray(profile.interests)) return [];
  return profile.interests;
}

// Onboarding chips + submit
(function initOnboarding(){
  const chipsWrap = document.getElementById("interestChips");
  const form = document.getElementById("onboardingForm");
  if (!chipsWrap || !form) return;

  const selected = new Set();

  chipsWrap.addEventListener("click", (e) => {
    const btn = e.target.closest(".chip");
    if (!btn) return;
    const val = btn.getAttribute("data-chip");
    if (!val) return;

    if (selected.has(val)) {
      selected.delete(val);
      btn.classList.remove("selected");
    } else {
      selected.add(val);
      btn.classList.add("selected");
    }
  });

  form.addEventListener("submit", (e) => {
    e.preventDefault();

    const age = document.getElementById("age").value.trim();
    const role = document.getElementById("role").value;
    const other = document.getElementById("interestOther").value.trim();
    const msg = document.getElementById("onboardingMsg");

    let interests = Array.from(selected);

    if (other) {
      const extra = other.split(",").map(s => s.trim()).filter(Boolean);
      interests = interests.concat(extra);
    }

    if (interests.length === 0) {
      if (msg) msg.textContent = "Pick at least 1 interest so we can personalize your feed.";
      return;
    }

    ssSetProfile({ age: Number(age), role, interests });
    window.location.href = "home.html";
  });
})();

// Home feed render
(function initHomeFeed(){
  const feedEl = document.getElementById("feed");
  const pillEl = document.getElementById("interestPill");
  if (!feedEl || !pillEl) return;

  const profile = ssGetProfile();
  if (!profile) {
    pillEl.textContent = "No profile found — tap Edit";
  } else {
    pillEl.textContent = `Interests: ${profile.interests.slice(0,3).join(", ")}${profile.interests.length > 3 ? " +" : ""}`;
  }

  const interests = ssGetSelectedInterests().map(s => s.toLowerCase());
  const papers = (window.SS_PAPERS || []).filter(p =>
    interests.length === 0 ? true : interests.includes(String(p.topic).toLowerCase())
  );

  if (papers.length === 0) {
    feedEl.innerHTML = `<p class="helper">No matching papers yet. Try adding more topics.</p>`;
    return;
  }

  feedEl.innerHTML = papers.map(p => `
    <div class="paperCard" role="button" tabindex="0" data-paper-id="${p.id}">
      <div class="paperTop">
        <div class="topicTag">${p.topic}</div>
        <div class="helper">${p.year}</div>
      </div>
      <div class="paperTitle">${p.title}</div>
      <p class="paperMeta">${p.authors} • ${p.journal}</p>
    </div>
  `).join("");

  feedEl.addEventListener("click", (e) => {
    const card = e.target.closest(".paperCard");
    if (!card) return;
    const id = card.getAttribute("data-paper-id");
    if (!id) return;
    localStorage.setItem("ss_selected_paper", id);
    window.location.href = "paper.html";
  });
})();

// Paper details page
(function initPaperPage(){
  const titleEl = document.getElementById("paperTitle");
  const metaEl = document.getElementById("paperMeta");
  const absEl = document.getElementById("paperAbstract");
  const openEl = document.getElementById("openPubMed");
  if (!titleEl || !metaEl || !absEl || !openEl) return;

  const id = localStorage.getItem("ss_selected_paper");
  const paper = (window.SS_PAPERS || []).find(p => p.id === id);

  if (!paper) {
    titleEl.textContent = "Paper not found";
    metaEl.textContent = "Go back and pick a paper from Home.";
    absEl.textContent = "";
    openEl.style.display = "none";
    return;
  }

  titleEl.textContent = paper.title;
  metaEl.textContent = `${paper.authors} • ${paper.journal} • ${paper.year} • ${paper.id}`;
  absEl.textContent = paper.abstract;
  openEl.href = paper.pubmedUrl || "#";
})();

// Settings reset
(function initSettings(){
  const btn = document.getElementById("resetProfile");
  const msg = document.getElementById("settingsMsg");
  if (!btn) return;

  btn.addEventListener("click", () => {
    localStorage.removeItem("ss_profile");
    localStorage.removeItem("ss_selected_paper");
    if (msg) msg.textContent = "Profile reset. Go to onboarding to set it up again.";
  });
})();
(function initLoginRedirect(){
  const form = document.getElementById("loginForm");
  if (!form) return;

  form.addEventListener("submit", (e) => {
    e.preventDefault();

    const msg = document.getElementById("loginMsg");
    if (msg) msg.textContent = "Logging in…";

    // Demo “success”: after real backend auth, replace this with your API call.
    setTimeout(() => {
      const profile = (() => {
        try { return JSON.parse(localStorage.getItem("ss_profile") || "null"); }
        catch { return null; }
      })();
        localStorage.setItem("ss_logged_in", "true");

      // ✅ After Log In → onboarding automatically (if profile missing)
      if (!profile) {
        window.location.href = "onboarding.html";
      } else {
        window.location.href = "home.html";
      }
    }, 250);
  });
})();
(function initSignupRedirect(){
  const form = document.getElementById("signupForm");
  if (!form) return;

  form.addEventListener("submit", (e) => {
    e.preventDefault();

    const msg = document.getElementById("signupMsg");
    if (msg) msg.textContent = "Creating account…";

    // Demo: backend will replace this
    setTimeout(() => {
      if (msg) msg.textContent = "✅ Account created. Redirecting to login…";
      setTimeout(() => (window.location.href = "login.html"), 350);
    }, 250);
  });
})();

(function authGuard(){
  const onHome = !!document.getElementById("feed");
  if (!onHome) return;

  const loggedIn = localStorage.getItem("ss_logged_in") === "true";
  if (!loggedIn) window.location.href = "login.html";
})();

function ssApplyPrefs(){
  const theme = localStorage.getItem("ss_theme") || "dark";
  document.documentElement.setAttribute("data-theme", theme);

  const font = localStorage.getItem("ss_font") || "system";
  const size = localStorage.getItem("ss_fontSize") || "16";

  const fontMap = {
    system: 'ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Arial',
    serif: 'ui-serif, Georgia, "Times New Roman", Times, serif',
    mono: 'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace'
  };

  document.documentElement.style.setProperty("--font-family", fontMap[font] || fontMap.system);
  document.documentElement.style.setProperty("--font-scale", `${size}px`);
}

ssApplyPrefs();

(function initSettingsPage(){
  const toggleTheme = document.getElementById("toggleTheme");
  const fontFamily = document.getElementById("fontFamily");
  const fontSize = document.getElementById("fontSize");
  const fontSizeValue = document.getElementById("fontSizeValue");
  const logoutBtn = document.getElementById("logoutBtn");
  const msg = document.getElementById("settingsMsg");

  if (!toggleTheme && !fontFamily && !fontSize && !logoutBtn) return;

  // set current values
  const currentFont = localStorage.getItem("ss_font") || "system";
  const currentSize = localStorage.getItem("ss_fontSize") || "16";
  if (fontFamily) fontFamily.value = currentFont;
  if (fontSize) fontSize.value = currentSize;
  if (fontSizeValue) fontSizeValue.textContent = `Current: ${currentSize}px`;

  if (toggleTheme) {
    toggleTheme.addEventListener("click", () => {
      const current = localStorage.getItem("ss_theme") || "dark";
      const next = current === "dark" ? "light" : "dark";
      localStorage.setItem("ss_theme", next);
      ssApplyPrefs();
      if (msg) msg.textContent = `Theme set to ${next}.`;
    });
  }

  if (fontFamily) {
    fontFamily.addEventListener("change", () => {
      localStorage.setItem("ss_font", fontFamily.value);
      ssApplyPrefs();
      if (msg) msg.textContent = "Font updated.";
    });
  }

  if (fontSize) {
    fontSize.addEventListener("input", () => {
      localStorage.setItem("ss_fontSize", fontSize.value);
      ssApplyPrefs();
      if (fontSizeValue) fontSizeValue.textContent = `Current: ${fontSize.value}px`;
    });
  }

  if (logoutBtn) {
    logoutBtn.addEventListener("click", () => {
      localStorage.removeItem("ss_logged_in");
      if (msg) msg.textContent = "Logged out.";
      setTimeout(() => (window.location.href = "login.html"), 250);
    });
  }
})();
(function initTopicsSelection(){
  const tiles = document.querySelectorAll(".topicTile");
  if (!tiles.length) return;

  const msg = document.getElementById("topicsMsg");

  // load current interests
  let profile = null;
  try { profile = JSON.parse(localStorage.getItem("ss_profile") || "null"); } catch {}
  if (!profile) profile = { interests: [] };
  if (!Array.isArray(profile.interests)) profile.interests = [];

  const interests = new Set(profile.interests);

  // mark selected
  tiles.forEach(t => {
    const topic = t.getAttribute("data-topic");
    if (topic && interests.has(topic)) t.classList.add("selected");
  });

  tiles.forEach(tile => {
    tile.addEventListener("click", () => {
      const topic = tile.getAttribute("data-topic");
      if (!topic) return;

      if (interests.has(topic)) {
        interests.delete(topic);
        tile.classList.remove("selected");
      } else {
        interests.add(topic);
        tile.classList.add("selected");
      }

      profile.interests = Array.from(interests);
      localStorage.setItem("ss_profile", JSON.stringify(profile));

      if (msg) msg.textContent = "Topics updated ✅ (Home feed can refresh based on these).";
    });
  });
})();
