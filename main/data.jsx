import { Event, trigger, everyTime, once } from "https://deno.land/x/good@0.5.15/events.js"

export const events = {
    annotationAdded: new Event([]),
    wigDataAdded: new Event([]),
}

export const data = {
    annotation: null,
    wigGroups: [],
    samples: [
        { disabled: false, condition: "Cholesterol", name: "cholesterol_H37Rv_rep1.wig", },
        { disabled: false, condition: "Cholesterol", name: "cholesterol_H37Rv_rep2.wig", },
        { disabled: false, condition: "Cholesterol", name: "cholesterol_H37Rv_rep3.wig", },
        { disabled: false, condition: "Glycerol",    name: "glycerol_H37Rv_rep1.wig",    },
        { disabled: false, condition: "Glycerol",    name: "glycerol_H37Rv_rep2.wig",    },
    ],
    conditions: [
        { disabled: false, name: "Cholesterol", },
        { disabled: false, name: "Cholesterol", },
        { disabled: false, name: "Cholesterol", },
        { disabled: false, name: "Glycerol",    },
        { disabled: false, name: "Glycerol",    },
    ],
    panelInfo: null,
}