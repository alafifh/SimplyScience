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
