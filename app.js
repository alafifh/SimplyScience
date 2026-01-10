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

