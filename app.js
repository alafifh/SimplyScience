const API_BASE = "https://simply-science.onrender.com/";

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

      // ✅ After Log In → onboarding automatically (if profile missing)
      if (!profile) {
        window.location.href = "onboarding.html";
      } else {
        window.location.href = "home.html";
      }
    }, 250);
  });
})();
async function fetchFeedFromBackend(query) {
  const res = await fetch(`${API_BASE}/feed`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      query: query,
      max_results: 5
    })
  });

  const data = await res.json();
  if (!data.ok) throw new Error("Feed failed");
  return data.facts;
}
document.getElementById("searchBtn").addEventListener("click", async () => {
  const q = document.getElementById("q").value.trim();
  if (!q) return;

  const facts = await fetchFeedFromBackend(q);

  // render facts into UI
  console.log(facts);
});
async function ssFetchFacts(query) {
  const res = await fetch(`${API_BASE}/feed`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ query, max_results: 5 })
  });

  // If CORS or server errors happen, this may throw
  const data = await res.json();
  if (!data.ok) {
    throw new Error(data.error || "Backend returned ok=false");
  }
  return data.facts || [];
}

(function initSearchPage() {
  const btn = document.getElementById("searchBtn");
  const input = document.getElementById("q");
  const status = document.getElementById("searchStatus");
  const results = document.getElementById("searchResults");

  if (!btn || !input || !results) return; // not on search page

  btn.addEventListener("click", async () => {
    const q = input.value.trim();
    results.innerHTML = "";
    if (!q) {
      if (status) status.textContent = "Type a topic (e.g., diabetes, brain cancer, sleep).";
      return;
    }

    try {
      if (status) status.textContent = "Searching PubMed + generating claims…";
      btn.disabled = true;

      const facts = await ssFetchFacts(q);

      if (!facts.length) {
        if (status) status.textContent = "No high-evidence claims returned. Try a broader query.";
        return;
      }

      if (status) status.textContent = `Results for: ${q}`;

      results.innerHTML = facts.map(f => `
        <div class="paperCard">
          <div class="paperTop">
            <div class="topicTag">${escapeHtml(f.category || "Claim")}</div>
            <div class="helper">${escapeHtml(f.evidence || "")}</div>
          </div>
          <div class="paperTitle">${escapeHtml(f.text || "")}</div>
          <p class="paperMeta">Claim ID: ${escapeHtml(f.id || "")}</p>
        </div>
      `).join("");
    } catch (err) {
      if (status) status.textContent = `Error: ${err.message}`;
    } finally {
      btn.disabled = false;
    }
  });

  // Enter key triggers search
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter") btn.click();
  });
})();

function escapeHtml(str) {
  return String(str)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}