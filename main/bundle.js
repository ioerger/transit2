// deno-fmt-ignore-file
// deno-lint-ignore-file
// This code was bundled using `deno bundle` and it's not recommended to edit it manually

var e = "", t = "";
function o(r1) {
    var p1, a1, l, s, c1 = arguments, i1 = this, n = 0, d = [], h = 0, u = [], f = 0;
    d.root = !0;
    var g = function(e1, o1, r2) {
        void 0 === o1 && (o1 = []);
        var p2 = 0;
        return (e1 = r2 || e1 !== t ? e1.replace(/\ue001/g, (e)=>u[f++]
        ) : u[f++].slice(1, -1)) ? (e1.replace(/\ue000/g, (t, r3)=>(r3 && o1.push(e1.slice(p2, r3)), p2 = r3 + 1, o1.push(c1[++h]))
        ), p2 < e1.length && o1.push(e1.slice(p2)), o1.length > 1 ? o1 : o1[0]) : e1;
    }, m = ()=>{
        [d, s, ...p1] = d, d.push(i1(s, ...p1));
    };
    return r1.join(e).replace(/<!--[^]*-->/g, "").replace(/<!\[CDATA\[[^]*\]\]>/g, "").replace(/('|")[^\1]*?\1/g, (e2)=>(u.push(e2), t)
    ).replace(/\s+/g, " ").replace(/(?:^|>)([^<]*)(?:$|<)/g, (e3, t1, r4, p3)=>{
        var c, i;
        if (r4 && p3.slice(n, r4).replace(/(\S)\/$/, "$1 /").split(" ").map((e4, t2)=>{
            if ("/" === e4[0]) c = i || e4.slice(1) || 1;
            else if (t2) {
                if (e4) {
                    var r5 = d[2] || (d[2] = {});
                    "..." === e4.slice(0, 3) ? Object.assign(r5, arguments[++h]) : ([a1, l] = e4.split("="), r5[g(a1)] = !l || g(l));
                }
            } else {
                for(i = g(e4); o.close[d[1] + i];)m();
                d = [
                    d,
                    i,
                    null
                ], o.empty[i] && (c = i);
            }
        }), c) for(m(); s !== c && o.close[s];)m();
        n = r4 + e3.length, t1 && " " !== t1 && g((s = 0, t1), d, !0);
    }), d.root || m(), d.length > 1 ? d : d[0];
}
o.empty = {}, o.close = {}, "area base br col command embed hr img input keygen link meta param source track wbr ! !doctype ? ?xml".split(" ").map((e5)=>o.empty[e5] = o.empty[e5.toUpperCase()] = !0
);
var r = {
    li: "",
    dt: "dd",
    dd: "dt",
    p: "address article aside blockquote details div dl fieldset figcaption figure footer form h1 h2 h3 h4 h5 h6 header hgroup hr main menu nav ol pre section table",
    rt: "rp",
    rp: "rt",
    optgroup: "",
    option: "optgroup",
    caption: "tbody thead tfoot tr colgroup",
    colgroup: "thead tbody tfoot tr caption",
    thead: "tbody tfoot caption",
    tbody: "tfoot caption",
    tfoot: "caption",
    tr: "tbody tfoot",
    td: "th tr",
    th: "td tr tbody"
}, p = function(e6) {
    [
        ...r[e6].split(" "),
        e6
    ].map((t3)=>{
        o.close[e6] = o.close[e6.toUpperCase()] = o.close[e6 + t3] = o.close[e6.toUpperCase() + t3] = o.close[e6 + t3.toUpperCase()] = o.close[e6.toUpperCase() + t3.toUpperCase()] = !0;
    });
};
for(var a in r)p(a);
const xhtm = o;
const kebabCase = (string)=>string.replace(/[a-z]([A-Z])(?=[a-z])/g, (each)=>`${each[0]}-${each.slice(1).toLowerCase()}`
    )
;
const isConstructor = (obj)=>!!obj.prototype && !!obj.prototype.constructor.name
;
const attachProperties = (source, target)=>{
    const attributes = Object.getOwnPropertyDescriptors(source);
    const propertiesDefition = {};
    for (const [key, value] of Object.entries(attributes)){
        if ([
            'constructor',
            'prototype',
            'length', 
        ].includes(key)) {
            continue;
        }
        propertiesDefition[key] = {
            get: ()=>ElementalClass[key]
        };
    }
    Object.defineProperties(target, propertiesDefition);
    return target;
};
class ElementalClass {
    constructor(components = {}, options = {}){
        const { middleware , errorComponentFactory  } = options || {};
        this.components = components || {};
        this.middleware = middleware || {};
        this.errorComponentFactory = errorComponentFactory || defaultErrorComponentFactory;
        this.html = this.createElement;
        this.xhtm = xhtm.bind((...args)=>this.createElement(...args)
        );
    }
    static debug = false;
    static allTags = Symbol.for("allTags");
    static exclusivelySvgElements = new Set([
        "svg",
        "animate",
        "animateMotion",
        "animateTransform",
        "circle",
        "clipPath",
        "defs",
        "desc",
        "discard",
        "ellipse",
        "feBlend",
        "feColorMatrix",
        "feComponentTransfer",
        "feComposite",
        "feConvolveMatrix",
        "feDiffuseLighting",
        "feDisplacementMap",
        "feDistantLight",
        "feDropShadow",
        "feFlood",
        "feFuncA",
        "feFuncB",
        "feFuncG",
        "feFuncR",
        "feGaussianBlur",
        "feImage",
        "feMerge",
        "feMergeNode",
        "feMorphology",
        "feOffset",
        "fePointLight",
        "feSpecularLighting",
        "feSpotLight",
        "feTile",
        "feTurbulence",
        "filter",
        "foreignObject",
        "g",
        "hatch",
        "hatchpath",
        "image",
        "line",
        "linearGradient",
        "marker",
        "mask",
        "mesh",
        "meshgradient",
        "meshpatch",
        "meshrow",
        "metadata",
        "mpath",
        "path",
        "pattern",
        "polygon",
        "polyline",
        "radialGradient",
        "rect",
        "set",
        "stop",
        "switch",
        "symbol",
        "text",
        "textPath",
        "tspan",
        "unknown",
        "use",
        "view", 
    ]);
    static randomId = (name)=>`${name}${Math.random()}`.replace(".", "")
    ;
    static appendChildren = function(element, ...children) {
        for (const each of children){
            if (typeof each == 'string') {
                element.appendChild(new window.Text(each));
            } else if (each == null) {
                element.appendChild(new window.Text(""));
            } else if (!(each instanceof Object)) {
                element.appendChild(new window.Text(`${each}`));
            } else if (each instanceof Node) {
                element.appendChild(each);
            } else if (each instanceof Array) {
                ElementalClass.appendChildren(element, ...each);
            } else if (each instanceof Function) {
                ElementalClass.appendChildren(element, each());
            } else if (each instanceof Promise) {
                const elementPromise = each;
                const placeholder = elementPromise.placeholder || document.createElement("div");
                setTimeout(async ()=>placeholder.replaceWith(await elementPromise)
                , 0);
                element.appendChild(placeholder);
            } else if (each != null && each instanceof Object) {
                element.appendChild(each);
            }
        }
        return element;
    };
    static css = function(first, ...args) {
        if (typeof first == 'string') {
            return first;
        } else if (first == null) {
            return "";
        } else if (first instanceof Array) {
            const strings = first;
            const values = args;
            let finalString = "";
            for (const each of strings){
                finalString += each;
                if (values.length > 0) {
                    if (value instanceof Object) {
                        finalString += Elemental.css(value);
                    } else {
                        finalString += `${values.shift()}`;
                    }
                }
            }
            return finalString;
        } else if (first instanceof Object) {
            let finalString = "";
            for (const [key, value] of Object.entries(first)){
                if (value != null) {
                    finalString += `${kebabCase(key)}: ${value};`;
                }
            }
            return finalString;
        } else {
            return first;
        }
    };
    static combineClasses = (...classes)=>{
        classes = classes.filter((each)=>each != null
        );
        let classesFinalList = [];
        for (let eachEntry of classes){
            if (typeof eachEntry == 'string') {
                eachEntry = eachEntry.split(" ");
            }
            if (eachEntry instanceof Array) {
                eachEntry = eachEntry.flat(Infinity);
                for (let eachName of eachEntry){
                    classesFinalList.push(eachName);
                }
            } else if (eachEntry instanceof Object) {
                for (const [className, enabled] of Object.entries(eachEntry)){
                    if (enabled) {
                        classesFinalList.push(className);
                    }
                }
            }
        }
        return classesFinalList;
    };
    createElement(...args) {
        if (args[0] instanceof Array) {
            return this.xhtm(...args);
        } else {
            ElementalClass.debug && console.debug(`args is:`, args);
            for (const middleware of (this.middleware[ElementalClass.allTags] || []).concat(this.middleware[args[0]] || [])){
                try {
                    args = eachMiddleWare(args);
                } catch (error) {
                    console.error("[ElementalClass] one of the middleware functions failed:", eachMiddleWare, args);
                }
            }
            let [key, properties, ...children] = args;
            ElementalClass.debug && console.debug(`key, properties, children is:`, key, properties, children);
            if (this.components[key] instanceof Function) {
                key = this.components[key];
            }
            if (key instanceof Function) {
                let output;
                try {
                    output = isConstructor(key) ? new key({
                        ...properties,
                        children
                    }) : key({
                        ...properties,
                        children
                    });
                } catch (error) {
                    return this.errorComponentFactory({
                        ...properties,
                        children
                    }, key, error);
                }
                if (output instanceof Promise) {
                    const elementPromise = output;
                    const placeholder = elementPromise.placeholder || document.createElement("div");
                    setTimeout(async ()=>placeholder.replaceWith(await elementPromise)
                    , 0);
                    return placeholder;
                } else {
                    return output;
                }
            }
            const isSvg = ElementalClass.exclusivelySvgElements.has(key);
            const element = isSvg ? document.createElementNS('http://www.w3.org/2000/svg', key) : document.createElement(key);
            if (properties instanceof Object) {
                for (let [key, value] of Object.entries(properties)){
                    if (key == 'class') {
                        if (value instanceof Array) {
                            value = value.join(" ");
                        } else if (value instanceof Object) {
                            let newValue = "";
                            for (const [classString, enable] of Object.entries(value)){
                                if (enable) {
                                    newValue += classString;
                                }
                            }
                            value = newValue;
                        }
                    }
                    if (key == 'style') {
                        value = ElementalClass.css(value);
                    }
                    try {
                        if (value instanceof Array) {
                            value = value.join(" ");
                        }
                        if (key.slice(0, 2) == 'on' && value instanceof Function) {
                            element.addEventListener(key.slice(2).toLowerCase(), value);
                        } else {
                            const attributeName = isSvg ? kebabCase(key) : key;
                            element.setAttribute(attributeName, value);
                        }
                    } catch (error) {
                        try {
                            element[key] = value;
                        } catch (error) {}
                    }
                }
            }
            return ElementalClass.appendChildren(element, ...children);
        }
    }
    extend(additionalComponents = {}, options = {}) {
        const { middleware , ...other } = options || {};
        return Elemental({
            ...this.components,
            ...additionalComponents
        }, {
            middleware: {
                ...this.middleware,
                ...middleware
            },
            ...other
        });
    }
}
const Elemental = (...args)=>{
    const elementalObject = new ElementalClass(...args);
    const createElementFunction = elementalObject.createElement.bind(elementalObject);
    attachProperties(ElementalClass, createElementFunction);
    attachProperties(Object.getPrototypeOf(elementalObject), createElementFunction);
    return createElementFunction;
};
attachProperties(ElementalClass, Elemental);
function defaultErrorComponentFactory({ children , ...properties }, key, error) {
    const element = document.createElement("div");
    const errorDetails = document.createElement("code");
    const childContainer = document.createElement("div");
    element.setAttribute('style', `
        all:              unset;
        display:          flex;
        flex-direction:   column;
        padding:          1.5rem;
        background-color: #f5a5a8;
        color:            white;
        font-family:      -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Oxygen-Sans,Ubuntu,Cantarell,"Helvetica Neue",sans-serif;
        font-size:        18px;
        font-weight:      400;
        overflow:         auto;
    `);
    element.innerHTML = `I'm sorry, there was an error when loading this part of the page 🙁 `;
    let errorElementPart;
    if (typeof key == 'string') {
        errorElementPart = `<${key} />`;
    } else {
        try {
            errorElementPart = `<${key.prototype.constructor.name} />`;
        } catch (error) {
            errorElementPart = `<${key} />`;
        }
    }
    let errorJsonObject = {};
    for (const [key1, value] of Object.entries(properties)){
        try {
            errorJsonObject[key1] = JSON.parse(JSON.stringify(value));
        } catch (error) {
            errorJsonObject[key1] = `${value}`;
        }
    }
    errorDetails.innerHTML = `tag: ${errorElementPart}\nproperties: ${JSON.stringify(errorJsonObject, 0, 4)}\nerror: ${error}`;
    errorDetails.setAttribute('style', `
        padding: 1rem;
        background-color: #161b22;
        color: #789896;
        white-space: pre;
        max-width: 85vw;
        overflow: auto;
    `);
    element.appendChild(errorDetails);
    childContainer.setAttribute('style', `
        all: unset
        display: flex
        flex-direction: column
        margin-top: 1.3rem
    `);
    ElementalClass.appendChildren(childContainer, children);
    element.appendChild(childContainer);
    return element;
}
try {
    const originalHead = document.head;
    Object.defineProperty(document, "head", {
        set: (element)=>ElementalClass.appendChildren(originalHead, ...element.childNodes)
        ,
        get: ()=>originalHead
        ,
        writable: true
    });
} catch (error) {}
ElementalClass.combineClasses;
const html = Elemental();
ElementalClass.css;
ElementalClass.allTags;
var e1 = "", t1 = "";
function o1(r1) {
    var p1, a1, l, s, c1 = arguments, i1 = this, n = 0, d = [], h = 0, u = [], f = 0;
    d.root = !0;
    var g = function(e11, o11, r2) {
        void 0 === o11 && (o11 = []);
        var p2 = 0;
        return (e11 = r2 || e11 !== t1 ? e11.replace(/\ue001/g, (e)=>u[f++]
        ) : u[f++].slice(1, -1)) ? (e11.replace(/\ue000/g, (t, r3)=>(r3 && o11.push(e11.slice(p2, r3)), p2 = r3 + 1, o11.push(c1[++h]))
        ), p2 < e11.length && o11.push(e11.slice(p2)), o11.length > 1 ? o11 : o11[0]) : e11;
    }, m = ()=>{
        [d, s, ...p1] = d, d.push(i1(s, ...p1));
    };
    return r1.join(e1).replace(/<!--[^]*-->/g, "").replace(/<!\[CDATA\[[^]*\]\]>/g, "").replace(/('|")[^\1]*?\1/g, (e2)=>(u.push(e2), t1)
    ).replace(/\s+/g, " ").replace(/(?:^|>)([^<]*)(?:$|<)/g, (e3, t11, r4, p3)=>{
        var c, i;
        if (r4 && p3.slice(n, r4).replace(/(\S)\/$/, "$1 /").split(" ").map((e4, t2)=>{
            if ("/" === e4[0]) c = i || e4.slice(1) || 1;
            else if (t2) {
                if (e4) {
                    var r5 = d[2] || (d[2] = {});
                    "..." === e4.slice(0, 3) ? Object.assign(r5, arguments[++h]) : ([a1, l] = e4.split("="), r5[g(a1)] = !l || g(l));
                }
            } else {
                for(i = g(e4); o1.close[d[1] + i];)m();
                d = [
                    d,
                    i,
                    null
                ], o1.empty[i] && (c = i);
            }
        }), c) for(m(); s !== c && o1.close[s];)m();
        n = r4 + e3.length, t11 && " " !== t11 && g((s = 0, t11), d, !0);
    }), d.root || m(), d.length > 1 ? d : d[0];
}
o1.empty = {}, o1.close = {}, "area base br col command embed hr img input keygen link meta param source track wbr ! !doctype ? ?xml".split(" ").map((e5)=>o1.empty[e5] = o1.empty[e5.toUpperCase()] = !0
);
var r1 = {
    li: "",
    dt: "dd",
    dd: "dt",
    p: "address article aside blockquote details div dl fieldset figcaption figure footer form h1 h2 h3 h4 h5 h6 header hgroup hr main menu nav ol pre section table",
    rt: "rp",
    rp: "rt",
    optgroup: "",
    option: "optgroup",
    caption: "tbody thead tfoot tr colgroup",
    colgroup: "thead tbody tfoot tr caption",
    thead: "tbody tfoot caption",
    tbody: "tfoot caption",
    tfoot: "caption",
    tr: "tbody tfoot",
    td: "th tr",
    th: "td tr tbody"
}, p1 = function(e6) {
    [
        ...r1[e6].split(" "),
        e6
    ].map((t3)=>{
        o1.close[e6] = o1.close[e6.toUpperCase()] = o1.close[e6 + t3] = o1.close[e6.toUpperCase() + t3] = o1.close[e6 + t3.toUpperCase()] = o1.close[e6.toUpperCase() + t3.toUpperCase()] = !0;
    });
};
for(var a1 in r1)p1(a1);
const xhtm1 = o1;
const kebabCase1 = (string)=>string.replace(/[a-z]([A-Z])(?=[a-z])/g, (each)=>`${each[0]}-${each.slice(1).toLowerCase()}`
    )
;
const isConstructor1 = (obj)=>!!obj.prototype && !!obj.prototype.constructor.name
;
const attachProperties1 = (source, target)=>{
    const attributes = Object.getOwnPropertyDescriptors(source);
    const propertiesDefition = {};
    for (const [key, value] of Object.entries(attributes)){
        if ([
            'constructor',
            'prototype',
            'length', 
        ].includes(key)) {
            continue;
        }
        propertiesDefition[key] = {
            get: ()=>ElementalClass1[key]
        };
    }
    Object.defineProperties(target, propertiesDefition);
    return target;
};
class ElementalClass1 {
    constructor(components = {}, options = {}){
        const { middleware , errorComponentFactory  } = options || {};
        this.components = components || {};
        this.middleware = middleware || {};
        this.errorComponentFactory = errorComponentFactory || defaultErrorComponentFactory1;
        this.html = this.createElement;
        this.xhtm = xhtm1.bind((...args)=>this.createElement(...args)
        );
    }
    static debug = false;
    static allTags = Symbol.for("allTags");
    static exclusivelySvgElements = new Set([
        "svg",
        "animate",
        "animateMotion",
        "animateTransform",
        "circle",
        "clipPath",
        "defs",
        "desc",
        "discard",
        "ellipse",
        "feBlend",
        "feColorMatrix",
        "feComponentTransfer",
        "feComposite",
        "feConvolveMatrix",
        "feDiffuseLighting",
        "feDisplacementMap",
        "feDistantLight",
        "feDropShadow",
        "feFlood",
        "feFuncA",
        "feFuncB",
        "feFuncG",
        "feFuncR",
        "feGaussianBlur",
        "feImage",
        "feMerge",
        "feMergeNode",
        "feMorphology",
        "feOffset",
        "fePointLight",
        "feSpecularLighting",
        "feSpotLight",
        "feTile",
        "feTurbulence",
        "filter",
        "foreignObject",
        "g",
        "hatch",
        "hatchpath",
        "image",
        "line",
        "linearGradient",
        "marker",
        "mask",
        "mesh",
        "meshgradient",
        "meshpatch",
        "meshrow",
        "metadata",
        "mpath",
        "path",
        "pattern",
        "polygon",
        "polyline",
        "radialGradient",
        "rect",
        "set",
        "stop",
        "switch",
        "symbol",
        "text",
        "textPath",
        "tspan",
        "unknown",
        "use",
        "view", 
    ]);
    static randomId = (name)=>`${name}${Math.random()}`.replace(".", "")
    ;
    static appendChildren = function(element, ...children) {
        for (const each of children){
            if (typeof each == 'string') {
                element.appendChild(new window.Text(each));
            } else if (each == null) {
                element.appendChild(new window.Text(""));
            } else if (!(each instanceof Object)) {
                element.appendChild(new window.Text(`${each}`));
            } else if (each instanceof Node) {
                element.appendChild(each);
            } else if (each instanceof Array) {
                ElementalClass1.appendChildren(element, ...each);
            } else if (each instanceof Function) {
                ElementalClass1.appendChildren(element, each());
            } else if (each instanceof Promise) {
                const elementPromise = each;
                const placeholder = elementPromise.placeholder || document.createElement("div");
                setTimeout(async ()=>placeholder.replaceWith(await elementPromise)
                , 0);
                element.appendChild(placeholder);
            } else if (each != null && each instanceof Object) {
                element.appendChild(each);
            }
        }
        return element;
    };
    static css = function(first, ...args) {
        if (typeof first == 'string') {
            return first;
        } else if (first == null) {
            return "";
        } else if (first instanceof Array) {
            const strings = first;
            const values = args;
            let finalString = "";
            for (const each of strings){
                finalString += each;
                if (values.length > 0) {
                    if (value instanceof Object) {
                        finalString += Elemental1.css(value);
                    } else {
                        finalString += `${values.shift()}`;
                    }
                }
            }
            return finalString;
        } else if (first instanceof Object) {
            let finalString = "";
            for (const [key, value] of Object.entries(first)){
                if (value != null) {
                    finalString += `${kebabCase1(key)}: ${value};`;
                }
            }
            return finalString;
        } else {
            return first;
        }
    };
    static combineClasses = (...classes)=>{
        classes = classes.filter((each)=>each != null
        );
        let classesFinalList = [];
        for (let eachEntry of classes){
            if (typeof eachEntry == 'string') {
                eachEntry = eachEntry.split(" ");
            }
            if (eachEntry instanceof Array) {
                eachEntry = eachEntry.flat(Infinity);
                for (let eachName of eachEntry){
                    classesFinalList.push(eachName);
                }
            } else if (eachEntry instanceof Object) {
                for (const [className, enabled] of Object.entries(eachEntry)){
                    if (enabled) {
                        classesFinalList.push(className);
                    }
                }
            }
        }
        return classesFinalList;
    };
    createElement(...args) {
        if (args[0] instanceof Array) {
            return this.xhtm(...args);
        } else {
            ElementalClass1.debug && console.debug(`args is:`, args);
            for (const middleware of (this.middleware[ElementalClass1.allTags] || []).concat(this.middleware[args[0]] || [])){
                try {
                    args = eachMiddleWare(args);
                } catch (error) {
                    console.error("[ElementalClass] one of the middleware functions failed:", eachMiddleWare, args);
                }
            }
            let [key, properties, ...children] = args;
            ElementalClass1.debug && console.debug(`key, properties, children is:`, key, properties, children);
            if (this.components[key] instanceof Function) {
                key = this.components[key];
            }
            if (key instanceof Function) {
                let output;
                try {
                    output = isConstructor1(key) ? new key({
                        ...properties,
                        children
                    }) : key({
                        ...properties,
                        children
                    });
                } catch (error1) {
                    return this.errorComponentFactory({
                        ...properties,
                        children
                    }, key, error1);
                }
                if (output instanceof Promise) {
                    const elementPromise = output;
                    const placeholder = elementPromise.placeholder || document.createElement("div");
                    setTimeout(async ()=>placeholder.replaceWith(await elementPromise)
                    , 0);
                    return placeholder;
                } else {
                    return output;
                }
            }
            const isSvg = ElementalClass1.exclusivelySvgElements.has(key);
            const element = isSvg ? document.createElementNS('http://www.w3.org/2000/svg', key) : document.createElement(key);
            if (properties instanceof Object) {
                for (let [key, value] of Object.entries(properties)){
                    if (key == 'class') {
                        if (value instanceof Array) {
                            value = value.join(" ");
                        } else if (value instanceof Object) {
                            let newValue = "";
                            for (const [classString, enable] of Object.entries(value)){
                                if (enable) {
                                    newValue += classString;
                                }
                            }
                            value = newValue;
                        }
                    }
                    if (key == 'style') {
                        value = ElementalClass1.css(value);
                    }
                    try {
                        if (value instanceof Array) {
                            value = value.join(" ");
                        }
                        if (key.slice(0, 2) == 'on' && value instanceof Function) {
                            element.addEventListener(key.slice(2).toLowerCase(), value);
                        } else {
                            const attributeName = isSvg ? kebabCase1(key) : key;
                            element.setAttribute(attributeName, value);
                        }
                    } catch (error) {
                        try {
                            element[key] = value;
                        } catch (error) {}
                    }
                }
            }
            return ElementalClass1.appendChildren(element, ...children);
        }
    }
    extend(additionalComponents = {}, options = {}) {
        const { middleware , ...other } = options || {};
        return Elemental1({
            ...this.components,
            ...additionalComponents
        }, {
            middleware: {
                ...this.middleware,
                ...middleware
            },
            ...other
        });
    }
}
const Elemental1 = (...args)=>{
    const elementalObject = new ElementalClass1(...args);
    const createElementFunction = elementalObject.createElement.bind(elementalObject);
    attachProperties1(ElementalClass1, createElementFunction);
    attachProperties1(Object.getPrototypeOf(elementalObject), createElementFunction);
    return createElementFunction;
};
attachProperties1(ElementalClass1, Elemental1);
function defaultErrorComponentFactory1({ children , ...properties }, key, error2) {
    const element = document.createElement("div");
    const errorDetails = document.createElement("code");
    const childContainer = document.createElement("div");
    element.setAttribute('style', `
        all:              unset;
        display:          flex;
        flex-direction:   column;
        padding:          1.5rem;
        background-color: #f5a5a8;
        color:            white;
        font-family:      -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Oxygen-Sans,Ubuntu,Cantarell,"Helvetica Neue",sans-serif;
        font-size:        18px;
        font-weight:      400;
        overflow:         auto;
    `);
    element.innerHTML = `I'm sorry, there was an error when loading this part of the page 🙁 `;
    let errorElementPart;
    if (typeof key == 'string') {
        errorElementPart = `<${key} />`;
    } else {
        try {
            errorElementPart = `<${key.prototype.constructor.name} />`;
        } catch (error) {
            errorElementPart = `<${key} />`;
        }
    }
    let errorJsonObject = {};
    for (const [key1, value] of Object.entries(properties)){
        try {
            errorJsonObject[key1] = JSON.parse(JSON.stringify(value));
        } catch (error) {
            errorJsonObject[key1] = `${value}`;
        }
    }
    errorDetails.innerHTML = `tag: ${errorElementPart}\nproperties: ${JSON.stringify(errorJsonObject, 0, 4)}\nerror: ${error2}`;
    errorDetails.setAttribute('style', `
        padding: 1rem;
        background-color: #161b22;
        color: #789896;
        white-space: pre;
        max-width: 85vw;
        overflow: auto;
    `);
    element.appendChild(errorDetails);
    childContainer.setAttribute('style', `
        all: unset
        display: flex
        flex-direction: column
        margin-top: 1.3rem
    `);
    ElementalClass1.appendChildren(childContainer, children);
    element.appendChild(childContainer);
    return element;
}
try {
    const originalHead = document.head;
    Object.defineProperty(document, "head", {
        set: (element)=>ElementalClass1.appendChildren(originalHead, ...element.childNodes)
        ,
        get: ()=>originalHead
        ,
        writable: true
    });
} catch (error3) {}
ElementalClass1.combineClasses;
const html1 = Elemental1();
const css = ElementalClass1.css;
ElementalClass1.allTags;
window.Elemental = Elemental1;
const randomId = Elemental1.randomId;
const combineClasses = Elemental1.combineClasses;
const translateAlignment = (name)=>{
    if (name == "top" || name == "left") {
        return "flex-start";
    } else if (name == "bottom" || name == "right") {
        return "flex-end";
    } else {
        return name;
    }
};
const classIds = {
    column: randomId(`column`),
    row: randomId(`row`),
    popUp: randomId(`popUp`),
    button: randomId(`button`),
    code: randomId(`code`)
};
document.body.appendChild(html1("style", null, `
    .${classIds.column} {
        display: flex;
        flex-direction: column;
        transition: all 0.5s ease-in-out 0s;
    }
    .${classIds.row} {
        display: flex;
        flex-direction: row;
        transition: all 0.5s ease-in-out 0s;
    }
    .${classIds.popUp} {
        align-items: center;
        justify-content: center;
        position: absolute;
        top:0;
        left: 0;
        z-index: 5000;
        display: flex;
        width: 100vw;
        height: 100vh;
    }
    .${classIds.button} {
        background: whitesmoke;
        position: relative;
        box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
        transition: all 0.5s ease-in-out 0s;
    }
    .${classIds.code} {
        white-space: pre;
        background: gray;
        position: relative;
        color: white;
        font-size: 11pt;
        padding: 0.5rem 0.7rem;
        border-radius: 0.4rem;
        transition: all 0.5s ease-in-out 0s;
    }
`));
const hoverStyleHelper = ({ element , hoverStyle  })=>{
    if (hoverStyle) {
        let hoverStyleAlreadyActive = false;
        const helper = document.createElement("div");
        const hoverStyleAsString = `${css(hoverStyle)}`;
        helper.style.cssText = hoverStyleAsString;
        const styleObject = {};
        const keys = Object.values(helper.style);
        for (const key1 of keys){
            styleObject[key1] = helper.style[key1];
        }
        const valuesBefore = {};
        element.addEventListener("mouseover", ()=>{
            if (!hoverStyleAlreadyActive) {
                hoverStyleAlreadyActive = true;
                for (const key of keys){
                    valuesBefore[key] = element.style[key];
                }
                element.style.cssText += hoverStyleAsString;
            }
        });
        element.addEventListener("mouseout", ()=>{
            if (hoverStyleAlreadyActive) {
                hoverStyleAlreadyActive = false;
                const style = element.style;
                const mixinStyleObject = {};
                for (const [key, value] of Object.entries(styleObject)){
                    if (style[key] == value) {
                        mixinStyleObject[key] = valuesBefore[key];
                    }
                }
                const mixinStyles = `${css(mixinStyleObject)}`;
                style.cssText += mixinStyles;
                Object.assign(style, mixinStyleObject);
            }
        });
    }
};
function Box({ children , style , hoverStyle , row , column =true , center =false , verticalAlignment =null , horizontalAlignment =null , onBlur , onChange , onClick , onContextMenu , onDblClick , onMouseDown , onMouseEnter , onMouseLeave , onMouseMove , onMouseOut , onMouseOver , onMouseUp , ...otherArgs }) {
    let justify, align, text, theClass;
    if (!row) {
        if (center) {
            verticalAlignment = verticalAlignment || "center";
            horizontalAlignment = horizontalAlignment || "center";
        } else {
            verticalAlignment = verticalAlignment || "top";
            horizontalAlignment = horizontalAlignment || "left";
        }
        justify = verticalAlignment;
        align = horizontalAlignment;
        text = horizontalAlignment;
        theClass = classIds.column;
    } else {
        if (center) {
            verticalAlignment = verticalAlignment || "center";
            horizontalAlignment = horizontalAlignment || "center";
        } else {
            verticalAlignment = verticalAlignment || "top";
            horizontalAlignment = horizontalAlignment || "left";
        }
        justify = horizontalAlignment;
        align = verticalAlignment;
        text = horizontalAlignment;
        theClass = classIds.row;
    }
    const element = html1("div", Object.assign({
        class: combineClasses(theClass, otherArgs.class),
        style: `justify-content: ${translateAlignment(justify)}; align-items: ${translateAlignment(align)}; text-align: ${text}; ${css(style)}; ${css(otherArgs)};`,
        onBlur: onBlur,
        onChange: onChange,
        onClick: onClick,
        onContextMenu: onContextMenu,
        onDblClick: onDblClick,
        onMouseDown: onMouseDown,
        onMouseEnter: onMouseEnter,
        onMouseLeave: onMouseLeave,
        onMouseMove: onMouseMove,
        onMouseOut: onMouseOut,
        onMouseOver: onMouseOver,
        onMouseUp: onMouseUp
    }, otherArgs), children);
    hoverStyleHelper({
        element,
        hoverStyle
    });
    return element;
}
const Column = Box;
const Row = (arg)=>Box({
        ...arg,
        row: true
    })
;
const Input = ({ children , style , checked , value , ...otherArgs })=>{
    const element = html1("input", Object.assign({
        class: otherArgs.class,
        style: `${css(style)}; ${css(otherArgs)};`
    }, otherArgs));
    Object.assign(element, {
        checked,
        value
    });
    return element;
};
const Code = ({ children , style , hoverStyle , onMouseOver , onMouseOut , onClick , ...otherArgs })=>{
    const element = html1("code", Object.assign({
        class: combineClasses(classIds.code, otherArgs.class),
        style: `${css(style)}; ${css(otherArgs)};`,
        onClick: onClick,
        onMouseOver: onMouseOver,
        onMouseOut: onMouseOut
    }, otherArgs), children);
    hoverStyleHelper({
        element,
        hoverStyle
    });
    return element;
};
const askForFiles = async ()=>{
    return new Promise((resolve, reject)=>{
        const cleanResolve = (returnValue)=>{
            resolve(returnValue);
            window.removeEventListener("focus", listener);
            document.body.removeChild(filePicker);
        };
        const listener = ()=>cleanResolve([])
        ;
        window.addEventListener("focus", listener);
        let filePicker = html1("input", {
            type: "file",
            onInput: (event)=>{
                cleanResolve(event.target.files);
            },
            onBlur: (event)=>{
                cleanResolve([]);
            },
            hidden: true
        });
        document.body.appendChild(filePicker);
        filePicker.click();
    });
};
function Button({ children , style , hoverStyle , row , column =true , center =false , verticalAlignment =null , horizontalAlignment =null , onBlur , onChange , onClick , onContextMenu , onDblClick , onMouseDown , onMouseEnter , onMouseLeave , onMouseMove , onMouseOut , onMouseOver , onMouseUp , ...otherArgs }) {
    let justify, align, text, theClass;
    if (!row) {
        if (center) {
            verticalAlignment = verticalAlignment || "center";
            horizontalAlignment = horizontalAlignment || "center";
        } else {
            verticalAlignment = verticalAlignment || "top";
            horizontalAlignment = horizontalAlignment || "left";
        }
        justify = verticalAlignment;
        align = horizontalAlignment;
        text = horizontalAlignment;
        theClass = classIds.column;
    } else {
        if (center) {
            verticalAlignment = verticalAlignment || "center";
            horizontalAlignment = horizontalAlignment || "center";
        } else {
            verticalAlignment = verticalAlignment || "top";
            horizontalAlignment = horizontalAlignment || "left";
        }
        justify = horizontalAlignment;
        align = verticalAlignment;
        text = horizontalAlignment;
        theClass = classIds.row;
    }
    const element = html1("button", Object.assign({
        class: combineClasses(theClass, otherArgs.class),
        style: `justify-content: ${translateAlignment(justify)}; align-items: ${translateAlignment(align)}; text-align: ${text}; ${css(style)}; ${css(otherArgs)};`,
        onBlur: onBlur,
        onChange: onChange,
        onClick: onClick,
        onContextMenu: onContextMenu,
        onDblClick: onDblClick,
        onMouseDown: onMouseDown,
        onMouseEnter: onMouseEnter,
        onMouseLeave: onMouseLeave,
        onMouseMove: onMouseMove,
        onMouseOut: onMouseOut,
        onMouseOver: onMouseOver,
        onMouseUp: onMouseUp
    }, otherArgs), children);
    hoverStyleHelper({
        element,
        hoverStyle
    });
    return element;
}
class Event extends Set {
}
const trigger = async (event, ...args)=>Promise.all([
        ...event
    ].map((each)=>each(...args)
    ))
;
const everyTime = (event)=>({
        then: (action)=>event.add(action)
    })
;
const EasyFilePicker = ({ children , onChange , defaultWidth ="16rem" , backgroundColor ="#939393" , ...otherArgs })=>{
    let codeElement;
    return html(Row, {
        style: `
            --default-width: ${defaultWidth};
            width: var(--default-width);
            max-width: var(--default-width);
            overflow-x: visible;
        `
    }, codeElement = html(Code, Object.assign({
        style: `
                    width: var(--default-width);
                    max-width: var(--default-width);
                    display: block;
                    overflow-x: auto;
                    background: ${backgroundColor};
                    position: relative;
                    box-shadow: 0 4px 5px 0 rgba(0,0,0,0), 0 1px 10px 0 rgba(0,0,0,0), 0 2px 4px -1px rgba(0,0,0,0);
                    margin-bottom: 0.5rem;
                `,
        hoverStyle: `
                    min-width: max-content;
                    overflow-x: visible;
                    box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
                `,
        onClick: async (event)=>{
            const files = await askForFiles();
            if (files && files.length > 0) {
                codeElement.innerText = files[0].name;
                onChange(files[0]);
            }
        }
    }, otherArgs), html(Row, {
        horizontalAlignment: "center",
        width: "100%"
    }, children)));
};
const events = {
    annotationAdded: new Event([]),
    wigDataAdded: new Event([])
};
const data = {
    annotation: null,
    wigGroups: [],
    samples: [
        {
            disabled: false,
            condition: "Cholesterol",
            name: "cholesterol_H37Rv_rep1.wig"
        },
        {
            disabled: false,
            condition: "Cholesterol",
            name: "cholesterol_H37Rv_rep2.wig"
        },
        {
            disabled: false,
            condition: "Cholesterol",
            name: "cholesterol_H37Rv_rep3.wig"
        },
        {
            disabled: false,
            condition: "Glycerol",
            name: "glycerol_H37Rv_rep1.wig"
        },
        {
            disabled: false,
            condition: "Glycerol",
            name: "glycerol_H37Rv_rep2.wig"
        }, 
    ],
    conditions: [
        {
            disabled: false,
            name: "Cholesterol"
        },
        {
            disabled: false,
            name: "Cholesterol"
        },
        {
            disabled: false,
            name: "Cholesterol"
        },
        {
            disabled: false,
            name: "Glycerol"
        },
        {
            disabled: false,
            name: "Glycerol"
        }, 
    ],
    panelInfo: null
};
const Annotation = ({ children , style ,  })=>{
    return html(Column, null, html("span", {
        class: "custom-header"
    }, "Annotation File"), html(Row, {
        padding: "1.2rem 1rem",
        horizontalAlignment: `space-between`,
        verticalAlignment: "center",
        background: "whitesmoke",
        border: "lightgray solid 1px",
        "min-width": "fit-content",
        "border-radius": "0.7rem"
    }, html(EasyFilePicker, {
        defaultWidth: "18rem",
        onChange: (file)=>{
            data.annotation = file.name;
            console.log(`data.annotation is:`, data.annotation);
            trigger(events.annotationAdded);
        }
    }, "[Click to add Annotation File]")));
};
const CombinedWigElement = ({ onLoaded  })=>{
    const data1 = {
        comWigFile: null,
        metadataFile: null
    };
    let onLoadedWasCalled = false;
    const checkData = (newData)=>{
        Object.assign(data1, newData);
        console.log(`wig data is:`, data1);
        if (!onLoadedWasCalled && data1.comWigFile && data1.metadataFile) {
            onLoadedWasCalled = true;
            onLoaded(data1);
        }
    };
    return html(Column, {
        style: `
            padding: 1rem;
            border: gray 1px solid;
            border-radius: 1rem;
            background: rgba(255,255,255,0.6);
            box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
            margin-bottom: 1rem;
        `
    }, html(EasyFilePicker, {
        defaultWidth: "17rem",
        onChange: (file)=>checkData({
                comWigFile: file
            })
    }, "[Click to add Combined Wig]"), html(EasyFilePicker, {
        defaultWidth: "17rem",
        onChange: (file)=>checkData({
                metadataFile: file
            })
    }, "[Click to add Metadata]"));
};
const WigLoader = ({ children , style ,  })=>{
    let wigLoader;
    const onLoaded = (newWigGroup)=>{
        data.wigGroups.push(newWigGroup);
        console.log(`events.wigDataAdded`);
        trigger(events.wigDataAdded);
        wigLoader.appendChild(html(CombinedWigElement, {
            onLoaded: onLoaded
        }));
    };
    return html(Column, {
        name: "WigLoader",
        height: "100%"
    }, html("span", {
        class: "custom-header"
    }, "Wig Files"), html(Column, {
        name: "WigLoader-Column",
        padding: "1.2rem 1rem",
        horizontalAlignment: `space-between`,
        "min-width": "21.5rem",
        "max-width": "21.5rem",
        background: "rgba(50, 134, 153, 0.2)",
        border: "lightgray solid 1px",
        "border-radius": "0.7rem",
        height: "100%",
        overflow: "visible"
    }, wigLoader = html(Column, {
        hoverStyle: `
                        min-width: 150vw;
                    `,
        overflow: "auto",
        height: "100%"
    }, html(CombinedWigElement, {
        onLoaded: onLoaded
    }))));
};
const standardColumnNames = [
    "disabled",
    "name",
    "condition", 
];
const createGridRowElements = ({ rowData , columns  })=>{
    console.debug(`rowData is:`, rowData);
    let rowElements = [
        html(Input, {
            type: "checkbox",
            checked: rowData.disabled,
            onChange: (event)=>{
                rowData.disabled = event.target.checked;
            }
        })
    ];
    for (const eachColumn of columns){
        if (eachColumn == "disabled") {
            continue;
        }
        const value = rowData[eachColumn];
        if (value != null && !(rowData[eachColumn] instanceof Object)) {
            rowElements.push(html("div", {
                class: "custom-grid-cell"
            }, " ", rowData[eachColumn], " "));
        } else {
            rowElements.push(html("div", {
                class: "custom-grid-cell"
            }, " "));
        }
    }
    return rowElements;
};
const GridElement = ({ columnNames , children  })=>{
    return html(Column, {
        padding: "1.2rem 1rem",
        background: "rgba(187, 225, 161, 0.28)",
        border: "lightgray solid 1px",
        "border-radius": "0.7rem",
        "flex-grow": "1",
        style: `
                display: grid;
                grid-template-columns: ${columnNames.map((each)=>"auto"
        ).join(" ")};
                width: 100%;
                grid-column-gap: 1px;
                grid-row-gap: 5px;
            `
    }, children);
};
const SampleFileTable = ({ children , style ,  })=>{
    let sampleFileElement, gridElement;
    everyTime(events.wigDataAdded).then(()=>{
        console.log(`everyTime(events.wigDataAdded)`);
        const columns = new Set(standardColumnNames);
        for (const eachSample of data.samples){
            for (const [key, value] of Object.entries(eachSample)){
                if (!(value instanceof Object)) {
                    columns.add(key);
                }
            }
        }
        const elements = [];
        const headerElements = [];
        for (const eachColumn of columns){
            headerElements.push(html(Row, {
                background: "gray",
                color: "white",
                padding: "0.4rem 0.8rem"
            }, eachColumn));
        }
        elements.push(headerElements);
        for (const eachSample1 of data.samples){
            elements.push(createGridRowElements({
                columns,
                rowData: eachSample1
            }));
        }
        const newGrid = html(GridElement, {
            columnNames: [
                ...columns
            ]
        }, elements.flat());
        gridElement.remove();
        sampleFileElement.appendChild(newGrid);
        gridElement = newGrid;
    });
    return sampleFileElement = html(Column, {
        "flex-grow": "1",
        width: "100%"
    }, html("span", {
        class: "custom-header"
    }, "Sample Files"), gridElement = html(GridElement, {
        columnNames: standardColumnNames
    }));
};
const standardColumnNames1 = [
    "disabled",
    "name", 
];
const createGridRowElements1 = ({ rowData , columns  })=>{
    console.debug(`rowData is:`, rowData);
    let rowElements = [
        html(Input, {
            type: "checkbox",
            checked: rowData.disabled,
            onChange: (event)=>{
                rowData.disabled = event.target.checked;
            }
        })
    ];
    for (const eachColumn of columns){
        if (eachColumn == "disabled") {
            continue;
        }
        const value = rowData[eachColumn];
        if (value != null && !(rowData[eachColumn] instanceof Object)) {
            rowElements.push(html("div", {
                class: "custom-grid-cell"
            }, " ", rowData[eachColumn], " "));
        } else {
            rowElements.push(html("div", {
                class: "custom-grid-cell"
            }, " "));
        }
    }
    return rowElements;
};
const GridElement1 = ({ columnNames , children  })=>{
    return html(Column, {
        padding: "1.2rem 1rem",
        background: "rgba(225, 224, 161, 0.28)",
        border: "lightgray solid 1px",
        "border-radius": "0.7rem",
        "flex-grow": "1",
        style: `
                display: grid;
                grid-template-columns: ${columnNames.map((each)=>"auto"
        ).join(" ")};
                width: 100%;
                grid-column-gap: 1px;
                grid-row-gap: 5px;
            `
    }, children);
};
const ConditionsTable = ({ children , style ,  })=>{
    let sampleFileElement, gridElement;
    everyTime(events.wigDataAdded).then(()=>{
        console.log(`everyTime(events.wigDataAdded)`);
        const columns = new Set(standardColumnNames1);
        for (const eachSample of data.conditions){
            for (const [key, value] of Object.entries(eachSample)){
                if (!(value instanceof Object)) {
                    columns.add(key);
                }
            }
        }
        const elements = [];
        const headerElements = [];
        for (const eachColumn of columns){
            headerElements.push(html(Row, {
                background: "gray",
                color: "white",
                padding: "0.4rem 0.8rem"
            }, eachColumn));
        }
        elements.push(headerElements);
        for (const eachCondition of data.conditions){
            elements.push(createGridRowElements1({
                columns,
                rowData: eachCondition
            }));
        }
        const newGrid = html(GridElement1, {
            columnNames: [
                ...columns
            ]
        }, elements.flat());
        gridElement.remove();
        sampleFileElement.appendChild(newGrid);
        gridElement = newGrid;
    });
    return sampleFileElement = html(Column, {
        "flex-grow": "1",
        width: "100%"
    }, html("span", {
        class: "custom-header"
    }, "Conditions"), gridElement = html(GridElement1, {
        columnNames: standardColumnNames1
    }));
};
const ParameterPanel = ({ children , style ,  })=>{
    return html(Column, {
        width: "100%",
        "flex-grow": "1"
    }, html("span", {
        class: "custom-header"
    }, "Analysis Parameters"), html(Column, {
        padding: "1.2rem 1rem",
        horizontalAlignment: `space-between`,
        verticalAlignment: "center",
        background: "whitesmoke",
        border: "lightgray solid 1px",
        "min-width": "fit-content",
        "border-radius": "0.7rem",
        width: "100%",
        "flex-grow": "1",
        gap: "2rem"
    }, html(Row, {
        width: "100%",
        horizontalAlignment: "space-between"
    }, "Pseudocount ", html(Input, {
        type: "number",
        value: "5"
    })), html(Row, {
        width: "100%",
        horizontalAlignment: "center"
    }, html(Button, {
        onClick: ()=>alert("Running ...")
        ,
        background: "rgb(202, 104, 104)",
        color: "white"
    }, "Run"))));
};
document.body.append(html(Row, {
    height: "100vh",
    width: "100vw",
    overflow: "hidden"
}, html("style", null, `
        .custom-header {
            font-weight: 500;
            font-size: 14pt;
            color: gray;
            margin-bottom: 2px;
            margin-left: 2px;
        }
    `), html(Column, {
    name: "FilePanel",
    width: "27rem",
    background: "white",
    padding: "3rem",
    height: "100%",
    "max-height": "100%"
}, html(Annotation, null), html(Row, {
    height: `3rem`
}), html(Row, {
    "flex-grow": "1",
    "min-height": "20rem",
    height: "100%"
}, html(WigLoader, null))), html(Column, {
    name: "DataPanel",
    width: "27rem",
    background: "white",
    padding: "3rem",
    height: "100%",
    "max-height": "100%",
    "flex-grow": "1"
}, html(SampleFileTable, null), html(Row, {
    height: `3rem`
}), html(ConditionsTable, null)), html(Column, {
    name: "ParameterPanel",
    width: "27rem",
    background: "white",
    padding: "3rem",
    height: "100%",
    "max-height": "100%"
}, html(ParameterPanel, null))));

