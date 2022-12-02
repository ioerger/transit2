let depthToShow = 3 // cant go below 2

const styleElement = document.createElement('style')
styleElement.innerHTML = `

body .wy-menu-vertical li.current a,
body .wy-menu-vertical li.current li.toctree-l2 a,
body .wy-menu-vertical li.current li.toctree-l3 a,
body .wy-menu-vertical li.current li.toctree-l4 a,
body .wy-menu-vertical li.current li.toctree-l5 a,
body .wy-menu-vertical li.current li.toctree-l6 a,
body .wy-menu-vertical li.current li.toctree-l7 a,
body .wy-menu-vertical li.current li.toctree-l8 a,
body .wy-menu-vertical li.current li.toctree-l9 a,
body .wy-menu-vertical li.current li.toctree-l10 a {
    color: rgb(0,0,0);
}


body .wy-menu-vertical .toctree-l2  a { color: rgba(255,255,255, 0.8); }
body .wy-menu-vertical .toctree-l3  a { color: rgba(255,255,255, 0.7); }
body .wy-menu-vertical .toctree-l4  a { color: rgba(255,255,255, 0.5); }
body .wy-menu-vertical .toctree-l5  a { color: rgba(255,255,255, 0.45); }
body .wy-menu-vertical .toctree-l6  a { color: rgba(255,255,255, 0.42); }
body .wy-menu-vertical .toctree-l7  a { color: rgba(255,255,255, 0.4); }
body .wy-menu-vertical .toctree-l8  a { color: rgba(255,255,255, 0.38); }
body .wy-menu-vertical .toctree-l9  a { color: rgba(255,255,255, 0.37); }
body .wy-menu-vertical .toctree-l10 a { color: rgba(255,255,255, 0.36); }

.wy-menu-vertical li a,
.wy-menu-vertical li.current a {
    opacity: 0.86;
    padding-left: 2rem;
}

.wy-menu-vertical li.toctree-l1 a,
.wy-menu-vertical li.toctree-l1.current > a {
    opacity: 0.86;
    padding-left: 2rem;
}

.wy-menu-vertical li.toctree-l2 a,
.wy-menu-vertical li.toctree-l2.current > a {
    opacity: 0.86;
    padding-left: 3rem;
}

.wy-menu-vertical li.toctree-l3 a,
.wy-menu-vertical li.toctree-l3.current a {
    opacity: 0.86;
    padding-left: 4rem;
}

.wy-menu-vertical li.toctree-l3 li a,
.wy-menu-vertical li.toctree-l3 li.current a {
    opacity: 0.86;
    padding-left: 5rem;
}

li.toctree-l1 > a {
    text-decoration: underline;
}

body li.toctree-l2 > a {
    color: whitesmoke;
}
`
depthToShow -= 1
while (--depthToShow) {
    styleElement.innerHTML += `
        body .wy-menu-vertical li.toctree-l${depthToShow} > ul {
            display: block !important;
        }
    `
}

document.head.appendChild(styleElement)

window.addEventListener('load', (_event) => {
    var menu = document.querySelector(".wy-menu ul li:first-child")
    recurse(menu)
});

/**
 * Given a Node, it recursively goes through every child and checks if the child is expandable, it
 * expands it unless it is already expanded.
 * 
 * @param {Node} node 
 */
function recurse(node) {
    if (is_expandable(node) && !is_expanded(node)) {
        node.classList.add("current")
    }

    // By default, children are not arrays, so we need to convert them
    children = Array.prototype.slice.call(node.children)

    children.forEach(recurse)
}

/**
 * Returns whether or not the given node is an expandable list.
 * 
 * @param {Node} node 
 * @returns {boolean} true if the node is a toctree that can be expanded, false otherwise.
 */
function is_expandable(node) {
    return node.className.includes("toctree-l")
}

/**
 * Returns whether or not the given expandable node is already expanded.
 * Nodes are considered expandaded if they are 'current'ly selected, so we take advantage of this.
 * 
 * @param {Node} node 
 * @returns {boolean} true if the node is already expanded, false otherwise.
 */
function is_expanded(node) {
    return node.classList.contains("current")
}