// Looping molecular movies (video.pz-loop) play muted and inline; a reader who asks for
// reduced motion gets the poster frame and play controls instead.
(function () {
    const query = window.matchMedia("(prefers-reduced-motion: reduce)");
    function apply() {
        document.querySelectorAll("video.pz-loop").forEach(function (video) {
            if (query.matches) {
                video.removeAttribute("autoplay");
                video.pause();
                video.controls = true;
            } else {
                video.controls = false;
                video.play().catch(function () {});
            }
        });
    }
    document.addEventListener("DOMContentLoaded", apply);
    query.addEventListener("change", apply);
})();
